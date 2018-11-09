#include<iostream>
#include<memory>
#include<vector>
#include<getopt.h>
#include<stdlib.h>
#include<set>
#include<string>
#include<math.h>

#include "white_alignment.h"

namespace white {

enum Direction { K_M, K_I, K_D, K_U }; //match, insertion, deletion, unknown

struct element {
	int score;
	Direction d;
};

struct before {
	int i;
	int j;
};

int max (int a, int b, int c)
{
	int rez;
	return c > (rez = a > b ? a : b) ? c : rez;
}

void max_MID (int match, int insertion, int deletion, element& e)
{
	int max = match;
	e.d = K_M;

	if (insertion > max)
	{
		max = insertion;
		e.d = K_I;
	}

	if (deletion > max)
	{
		max = deletion;
		e.d = K_D;
	}

	e.score = max;
}

void initialization_global (element* matrix, int query_length, int target_length, int gap)
{

	matrix[0].score = 0;
        matrix[0].d = K_U;

	for (int i = 1; i < target_length + 1; i++)
	{
               	matrix[i].score = gap * i;
		matrix[i].d = K_I;
	}

	for (int i = 1; i < query_length + 1; i++)
	{
		matrix[i * (target_length + 1)].score = gap * i;
		matrix[i * (target_length + 1)].d = K_D;
	}
}

void initialization_local_semi_global (element* matrix, int query_length, int target_length)
{

		matrix[0].score = 0;
		matrix[0].d = K_U;

		for (int i = 1; i < target_length + 1; i++)
                {
                	matrix[i].score = 0;
                	matrix[i].d = K_I;
        	}

                for (int i = 1; i < query_length + 1; i++)
               {
                        matrix[i * (target_length + 1)].score = 0;
                        matrix[i * (target_length + 1)].d = K_D;
               }
}

void print_matrix (element* matrix, int rows, int columns)
{
	for (int i = 0; i < rows; i++)
	{
        	for (int j = 0; j < columns; j++)
			printf ("%d_%d\t", matrix[i * columns + j].score, matrix[i * columns + j].d);

		printf ("\n");
	}
}

void matrix_calculation (	const char* query, unsigned int query_length,
				const char* target, unsigned int target_length,
				int match,
				int mismatch,
				int gap,
				element* matrix)
{
	for (unsigned int i = 1; i < query_length + 1; i++)
                for (unsigned int j = 1; j < target_length + 1; j++)
                {
                        int  mm = query[i-1] == target[j-1] ? match : mismatch;

                        int match = matrix[(i-1) * (target_length + 1) + (j-1)].score + mm;
                        int insertion = matrix[(i) * (target_length + 1) + (j-1)].score + gap;
                        int deletion = matrix[(i-1) * (target_length + 1) + (j)].score + gap;

                        max_MID (match, insertion, deletion, matrix[i * (target_length + 1) + j]);
                }

//	print_matrix (matrix, query_length + 1, target_length + 1);

}

void matrix_calculation_local ( const char* query, unsigned int query_length,
                                const char* target, unsigned int target_length,
                                int match,
                                int mismatch,
                                int gap,
				int& row,
				int& column,
                                element* matrix)
{
	int max = 0;

	 for (unsigned int i = 1; i < query_length + 1; i++)
                for (unsigned int j = 1; j < target_length + 1; j++)
                {
                        int  mm = query[i-1] == target[j-1] ? match : mismatch;

                        int match = matrix[(i-1) * (target_length + 1) + (j-1)].score + mm;
                        int insertion = matrix[(i) * (target_length + 1) + (j-1)].score + gap;
                        int deletion = matrix[(i-1) * (target_length + 1)+ (j)].score + gap;

			if (match < 0 && insertion < 0 && deletion < 0)
        		{
                		matrix[i * (target_length + 1) + j].score = 0;
                		matrix[i * (target_length + 1) + j].d = K_U;
        		}

        		else max_MID (match, insertion, deletion, matrix[i * (target_length + 1) + j]);

			if (max < matrix[i * (target_length + 1)+ j].score)
			{
				max = matrix[i * (target_length + 1)+ j].score;
				row = i;
				column = j;
			}
                }
//	print_matrix (matrix, query_length + 1, target_length + 1);
}

void write_int_into_char_array_reversed (int number, char* array, int& starting_position)
{

	                        std::string temp_number = std::to_string(number);
                                int size = temp_number.size();

                                const char* temp = temp_number.c_str();

                                for (int i = 0; i < size; i++)
                                        array[starting_position++] = temp[size - 1 - i];

}

void cigar_string (int query_length, int target_length, AlignmentType type, element* matrix, std::string& cigar, unsigned int& target_begin, int row, int column)
{
	char* r_cigar = new char[(query_length + target_length) * 2 + 1]; //cigar će biti generiran u obrnutom redoslijedu (reverse_cigar)

	int k = 0; //brojač koliko r_cigar ima elemenata
        char last = ' '; //zadnji element koji se pojavio
        int counter = 0; //brojač koliko je bilo uzastopce istih elemenata
        target_begin = 1; //defaultno 1

	int i = row;
        int j = column;
	before b;

	switch (type) {
		case global:
			break;

		case semi_global:
			if (column == target_length && row < query_length)
			{
				r_cigar[k++] = 'S';
				write_int_into_char_array_reversed (query_length - row, r_cigar, k);
			}
			break;

		case local:
			if (row < query_length)
			{
				r_cigar[k++] = 'S';
				write_int_into_char_array_reversed (query_length - row, r_cigar, k);
			}
			break;
	}

	while (1)
        {
                if (matrix[i * (target_length + 1) + j].d == K_M)
                {
                        if (last == ' ')
			{
				last = 'M';
	                 	counter++;
                        }

			else if (last == 'M')
			{
				counter++;
			}
			else
                        {
				r_cigar[k++] = last;

				write_int_into_char_array_reversed (counter, r_cigar, k);

				last = 'M';
                                counter = 1;
                        }

			b.i = i;
			b.j = j;
                        i--;
                        j--;

		}

		else if (matrix[i * (target_length + 1) + j].d == K_I)
                {
			if (last == ' ')
                        {
                                last = 'I';
                                counter++;
                        }

                        else if (last == 'I')
                        {
                                counter++;
                        }
                        else
                        {
                                r_cigar[k++] = last;

				write_int_into_char_array_reversed (counter, r_cigar, k);

                                last = 'I';
                                counter = 1;
                        }

			b.i = i;
                        b.j = j;
                        j--;

                }

		else if(matrix[i * (target_length + 1) + j].d == K_D)
		{
			if (last == ' ')
                        {
                                last = 'D';
                                counter++;
                        }

                        else if (last == 'D')
                        {
                                counter++;
                        }
                        else
                        {
                                r_cigar[k++] = last;

				write_int_into_char_array_reversed (counter, r_cigar, k);

                                last = 'D';
                                counter = 1;
                        }

			b.i = i;
                        b.j = j;
                        i--;
                }

		if (type == global)
		{
			if (i == 0 && j == 0)
			{
				r_cigar[k++] = last;
                        	write_int_into_char_array_reversed (counter, r_cigar, k);
				break;
			}
		}

		else if (type == semi_global)
		{
			if (i == 0)
			{
                                r_cigar[k++] = last;
				write_int_into_char_array_reversed (counter, r_cigar, k);

				target_begin = b.j;
				break;
			}
			else if (j == 0)
			{
                                r_cigar[k++] = last;
				write_int_into_char_array_reversed (counter, r_cigar, k);

				if (b.i - 1 > 0)
				{
					r_cigar[k++] = 'S';
					write_int_into_char_array_reversed (b.i - 1, r_cigar, k);
				}

                                target_begin = 1;
				break;
			}
		}

		else if (type == local)
		{
			if (matrix[i * (target_length + 1) + j].score == 0)
			{
                                r_cigar[k++] = last;
				write_int_into_char_array_reversed (counter, r_cigar, k);

				if (b.i - 1 > 0)
				{
                                	r_cigar[k++] = 'S';
					write_int_into_char_array_reversed (b.i - 1, r_cigar, k);
				}
				target_begin = b.j;
                                break;
			}
		}

	}

        //obrnimo dobiveni r_cigar
	r_cigar[k] = '\0';
        char* cigar_char = new char[k + 1];

        for (int i = 0; i < k; i++)
                cigar_char[i] = r_cigar[k-1-i];

	cigar_char[k] = '\0';
	cigar = std::string (cigar_char);

	delete[] cigar_char;
	delete[] r_cigar;
}

int global_alignment (	const char* query, unsigned int query_length,
			const char* target, unsigned int target_length,
			int match,
			int mismatch,
			int gap,
			std::string& cigar,
			unsigned int& target_begin)
{

	element* matrix = new element[(query_length + 1) * (target_length + 1)];

	initialization_global (matrix, query_length, target_length, gap);

        matrix_calculation (query, query_length, target, target_length, match, mismatch, gap, matrix);

	cigar_string (query_length, target_length, white::global, matrix, cigar, target_begin, query_length, target_length);

	int result = matrix[query_length * (target_length + 1) + target_length].score;

	delete[] matrix;

	return result;
}

int global_alignment (const char* query, unsigned int query_length,
                      const char* target, unsigned int target_length,
                      int match,
                      int mismatch,
                      int gap)
{
        std::string cigar;
        unsigned int target_begin;

        return global_alignment (query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);
}

int max_matrix_LR_LC (element* matrix, unsigned int query_length, unsigned int target_length, int& row, int& column)
{

	int max = matrix[query_length * (target_length + 1)].score;

	for (unsigned int j = 0; j < target_length + 1; j++)
		if (matrix[query_length * (target_length + 1) + j].score > max)
		{
			max = matrix[query_length * (target_length + 1) + j].score;
			row = query_length;
			column = j;
		}

	 for (unsigned int i = 0; i < query_length; i++)
                if (matrix[i * (target_length + 1) + target_length].score > max)
                {
                        max = matrix[i * (target_length + 1) + target_length].score;
                        row = i;
                        column = target_length;
                }

	return max;

}

int semi_global_alignment (     const char* query, unsigned int query_length,
                                const char* target, unsigned int target_length,
                                int match,
                                int mismatch,
                                int gap,
				std::string& cigar,
				unsigned int& target_begin)
{
	element* matrix = new element[(query_length + 1) * (target_length + 1)];

	int row = 0;
	int column = 0;

	initialization_local_semi_global (matrix, query_length, target_length);

        matrix_calculation (query, query_length, target, target_length, match, mismatch, gap, matrix);

	max_matrix_LR_LC (matrix, query_length, target_length, row, column);

	cigar_string (query_length, target_length, white::semi_global, matrix, cigar, target_begin, row, column);

	int result = matrix[row * (target_length + 1) + column].score;

	delete[] matrix;

	return result;
}

int semi_global_alignment (     const char* query, unsigned int query_length,
                                const char* target, unsigned int target_length,
                                int match,
                                int mismatch,
                                int gap)
{
        std::string cigar;
        unsigned int target_begin;

        return semi_global_alignment (query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);
}

int local_alignment (   const char* query, unsigned int query_length,
                        const char* target, unsigned int target_length,
                        int match,
                        int mismatch,
                        int gap,
			std::string& cigar,
                        unsigned int& target_begin)
{
	element* matrix = new element[(query_length + 1) * (target_length + 1)];

        int row = 0;
        int column = 0;

        initialization_local_semi_global (matrix, query_length, target_length);

        matrix_calculation_local (query, query_length, target, target_length, match, mismatch, gap, row, column, matrix);

	cigar_string (query_length, target_length, white::local, matrix, cigar, target_begin, row, column);

	int result = matrix[row*(target_length + 1) + column].score;

	delete[] matrix;

        return result;
}

int local_alignment (   const char* query, unsigned int query_length,
                        const char* target, unsigned int target_length,
                        int match,
                        int mismatch,
                        int gap)
{
        std::string cigar;
        unsigned int target_begin;

        return local_alignment (query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);
}

//glavna funkcija

int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap)
{
	std::string cigar;
        unsigned int target_begin;

	return pairwise_alignment (query, query_length, target, target_length, type, match, mismatch, gap, cigar, target_begin);
}

//nadjačana glavna funkcija

int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap,
                       std::string& cigar,
                       unsigned int& target_begin)
{
	 switch (type)
        {
                case global:
                        return global_alignment (query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);

                case semi_global:
                        return semi_global_alignment (query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);

                case local:
                        return local_alignment (query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);
		default :
			return 0;
        }

}
}
