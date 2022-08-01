


//   The C++ version of the Math, Number, Matrix, and Convert classes
//
//   Every line of the C++ code has to be compared to the Java code because
//   some errors were corrected in the Java version but not in the C++ version




#include <string>
#include <vector>
#include <iostream>
#include <cctype>
#include <thread>
#include <mutex>

// #include <cmath>
// #include <algorithm>


using namespace std;





class Number
{


	public:
	
	       ~Number();
		Number();
		
		Number(int i);
		Number(long l);
		Number(double d);
		Number(string digits);
		Number(string digits, int radix);
		Number(vector<int> vec_int);
		Number(vector<int> vec_int, bool signed1);
		
		Number(int real, int imag);
		Number(double real, double imag);
		Number(Number& real, Number& imag);
		
		Number(Number& original); // copy constructor
	
	
	private:
	
		static Number* double_to_number(double d);
		
		static Number* string_to_number( string &digits, int radix);
		static Number* string_to_number1(string &digits, int radix);
		
		Number& set_result(Number& number);
	
	
	public:
	
		Number& abs();
		
		Number& add(int addend);
		Number& add(long addend);
		Number& add(double addend);
		Number& add(Number& addend);
		
		Number& add_bit(long bit);
		
		Number& and_(Number& number);
		
		Number& arc_cos();
		Number& arc_sin();
		Number& arc_tan();
		
		Number& arc_cosh();
		Number& arc_sinh();
		Number& arc_tanh();
		
		long bit_count();
		
		Number& clear_bit(long bit);
		
		int compare(double val);
		
		int compare(int val);
		
		int compare(long val);
		
		int compare(Number& val);
		
		Number& complement();
		
		Number& cos();
	
	
	
	private:
	
		Number *cos_sin();
	
	
	public:
	
		Number& cosh();
		
		Number& cube();
		
		Number& divide(int divisor);
		Number& divide(long divisor);
		Number& divide(double divisor);
		Number& divide(Number& divisor);
	
	
	private:
	
		Number& divide_by_int(int divisor);
	
	
	public:
	
		double double_value();
		
		static Number& e_();
		static Number& e_(int digits);
		
		bool equals(int val);
		bool equals(long val);
		bool equals(double val);
		bool equals(Number& val);
		
		static Number& exp(int x);
		static Number& exp(Number& x);
		
		Number& factor();
		
		Number& factorial();
		
		Number& flip_bit(long bit);
		
		float float_value();
		
		Number& gcd(long n);
		
		Number& gcd(Number& n);
		
		int get_bit(long bit);
		
		long get_lowest_set_bit();
		
		int get_precision();
		
		int int_value();
		
		Number& inverse();
		
		static bool is_base_16(string &str);
		static bool is_base_64(string &str);
		
		bool is_complex();
		
		bool is_coprime_with(long n);
		bool is_coprime_with(Number& n);
		
		static bool is_digit_string(string &str, int radix);
		
		bool is_divisible_by(int n);
		bool is_divisible_by(Number& n);
		
		bool is_even();
		
		bool is_factorable(int max_prime);
		
		bool greater_than(int val);
		bool greater_than(double val);
		bool greater_than(long val);
		bool greater_than(Number& val);
		
		bool is_integer();
		
		static bool is_integer_string(string& str, int radix);
		
		bool less_than(double val);
		bool less_than(int val);
		bool less_than(long val);
		bool less_than(Number& val);
		
		static bool is_number_string(string& str, int radix);
		
		bool is_generator(Number& p);
		
		bool is_power(int n);
		
		bool is_power_of(int n);
		
		bool is_prime();
		bool is_probable_prime();
		bool is_probable_prime(int certainty);
		
		bool is_quadratic_residue(Number& p);
		
		bool is_square();
		
		Number& lambda();
		
		Number& lcm(long n);
		Number& lcm(Number& a);
		
		int length();
		
		Number& ln();
		
		long log2();
		
		Number& log();
		Number& log(int base);
		Number& log(Number& base);
	
	
	
	private:
	
		Number& log1();
	
	
	public:
	
		long long_value();
		
		Number& lvalue();
		
		Number& max(Number& val);
		
		Number& min(Number& val);
		
		Number& mod(double val);
		Number& mod(int val);
		Number& mod(long val);
		Number& mod(Number& mod);
		Number& mod(Number& n, Number& inv);
		
		Number& mod_divide(int divisor, Number& n);
		Number& mod_divide(Number& divisor, Number& n);
		
		Number& mod_inverse(long modulus);
		Number& mod_inverse(Number& n);
	
	
	
	private:
	
		Number& mod_inverse1(Number& n);
		
		Number& mod_multiply(Number& ar, Number& n, Number& n1, Number& r);
		
		Number& mod_pow1(Number& exp, Number& n);
	
	
	public:
	
		Number& mod_pow(int exp, int m);
		Number& mod_pow(int exp, Number& m);
		Number& mod_pow(Number& exp, int m);
		Number& mod_pow(Number& exp, Number& n);
		
		Number& mod_root(int k, Number& n);
		
		Number& multiply(int multiplier);
		Number& multiply(long multiplier);
		Number& multiply(double multiplier);
		Number& multiply(Number& multiplier);
		
		Number& negate();
		Number& negate(int n);
		Number& negate(Number& n);
		
		Number& next_prime();
		
		static int parse_int(string & s);
		static int parse_int(string & s, int radix);
		
		Number& phi();
		
		static Number& pi();
		static Number& pi(int digits);
		static Number& pi1(int digits);
		static Number& pi11(int digits);
		
		Number& pow(int exp);
		Number& pow(double exp);
		Number& pow(Number& exp);
	
	
	private:
	
		Number& pow1(Number& exp);
		
		Number& quad_divide(Number& divisor);
	
	
	public:
	
		// static Number& random(int digits, int radix);
		
		Number& root(int k);
		
		Number& round();
		
		Number& set_bit(long bit);
		
		Number& set_precision(int precision);
		
		Number& shift_left(long bits);
		Number& shift_left(long expansion, long bits);
		Number& shift_right(long bits);
		
		int signum();
		
		Number& sin();
		
		Number& sinh();
		
		Number& sqrt();
		
		Number& square();
		
		Number& subtract(int subtrahend);
		Number& subtract(long subtrahend);
		Number& subtract(double subtrahend);
		Number& subtract(Number& subtrahend);
		
		Number& tan();
		
		Number& tanh();
		
		bool test_bit(long bit);
		
		string to_alphabetical_string();
		
		Number& to_fraction();
		
		Number& to_imag();
		
		vector<int> to_int_array();
		vector<int> to_int_array(int digits, long radix);
		vector<int> to_int_vector();
		vector<int> to_int_vector(int digits);
		
		Number& to_integer();
		
		Number& to_real();
		
		vector<int> to_signed_int_vector(int length);
		
		string to_string();
		string to_string(int radix);
		string to_string(int digits, int radix);
	
	
	private:
	
		string to_string1(int radix);
		string to_string2(int length, int radix);
	
	
	public:
	
		Number& trim();
		Number& trim(int bits);
	
	
	
	private:
	
		vector<Number*> result;
		
		vector<int> vec_int = vector<int>(0);
		
		int intpoint = 0;
		
		int precision = 0;
		
		char sign = '+';
		
		
		
		// imaginary number
		
		vector<int> vec_int1 = vector<int>(0);;
		
		int intpoint1 = 0;
		
		int precision1 = 0;
		
		char sign1 = '+';
	
	
	
	
	//  Define the overloaded assignment operator
	
	public:
	
		Number& operator = (Number &original);
	
	
	
	
	
	//  Define the overloaded operators
	//
	//  +   -   *   /   ^   <<   >>   ...
	//
	//  +=  -=  *=  /=  ^=  <<=  >>=  ...
	
	
	public:
	
	
		Number& operator + (int b)
		{
			return this->add(b);
		}
		
		Number& operator - (int b)
		{
			return this->subtract(b);
		}
		
		Number& operator * (int b)
		{
			return this->multiply(b);
		}
		
		Number& operator / (int b)
		{
			return this->divide(b);
		}
		
		
		
		Number& operator + (Number& b)
		{
			return this->add(b);
		}
		
		Number& operator - (Number& b)
		{
			return this->subtract(b);
		}
		
		Number& operator * (Number& b)
		{
			return this->multiply(b);
		}
		
		Number& operator / (Number& b)
		{
			return this->divide(b);
		}
		
		
		
		
		Number& operator += (Number& b)
		{
			Number *sum = & this->add(b);
			
			this->vec_int = sum->vec_int;
			
			this->intpoint  = sum->intpoint;
			this->precision = sum->precision;
			this->sign      = sum->sign;
			
			delete sum;
			
			return *this;
		}
		
		
		Number& operator -= (Number& b)
		{
			Number *diff = & this->subtract(b);
			
			this->vec_int = diff->vec_int;
			
			this->intpoint  = diff->intpoint;
			this->precision = diff->precision;
			this->sign      = diff->sign;
			
			delete diff;
			
			return *this;
		}
		
		
		Number& operator *= (Number& b)
		{
			Number *product = & this->multiply(b);
			
			this->vec_int = product->vec_int;
			
			this->intpoint  = product->intpoint;
			this->precision = product->precision;
			this->sign      = product->sign;
			
			delete product;
			
			return *this;
		}
		
		
		Number& operator /= (Number& b)
		{
			Number *quotient = & this->divide(b);
			
			this->vec_int = quotient->vec_int;
			
			this->intpoint  = quotient->intpoint;
			this->precision = quotient->precision;
			this->sign      = quotient->sign;
			
			delete quotient;
			
			return *this;
		}
		
		
		//	...
		
		//	...
};


//  End class Number















class Matrix
{


	public:
	
	       ~Matrix();
		Matrix();
		
		Matrix(int rows, int columns);
		
		Matrix(vector<int> matrix, int rows, int columns);
		Matrix(vector<double> matrix, int rows, int columns);
		Matrix(vector<Number*> matrix, int rows, int columns);
		
		Matrix(vector<int> row);
		Matrix(vector<double> row);
		Matrix(vector<Number*> row);
		
		Matrix(string &str, int radix, int rows, int columns);
		
		Matrix(Matrix& original);  // copy constructor
		
	private:
	
		Matrix& set_result(Matrix& matrix);
	
	
	public:
	
		Matrix& abs();
		
		Matrix& add(int number);
		Matrix& add(Number& number);
		Matrix& add(Matrix& matrix);
		
		Matrix& and_(Number n);
		
		Matrix& augment(vector<int> b);
		Matrix& augment(vector<Number*> B);
		Matrix& augment(Matrix& B);
		
		int column_count();
		
		Matrix& delete_column(int column);
		
		Matrix& delete_row(int row);
		
		Number& determinant();
		Number& determinant(Number& modulus);
		Number& determinant1();
		Number& determinant1(Number& modulus);
		
		Matrix& divide(int n);
		Matrix& divide(Number& number);
		
		bool equals(Matrix& matrix);
		
		Number& get(int i, int j);
		Matrix& get(int row, int column, int rows, int columns);
		
		vector<Number*> get_column(int column);
		
		vector<Number*> get_diagonals();
		
		int get_precision();
		
		vector<Number*> get_row(int row);
		
		Matrix& identity_matrix(int size);
		
		Matrix& inverse();
		
		bool is_complex();
		
		bool is_echelon_form();
		
		bool is_identity_matrix();
		
		bool is_row_canonical_form();
		
		bool is_singular();
		bool is_singular(Number& n);
		
		bool is_square();
		
		Matrix *LU();
		Matrix *LU(Number& modulus);
		
		
		Matrix& mod(int modulus);
		Matrix& mod(Number& modulus);
		
		Matrix& mod_divide(int divisor, Number& n);
		Matrix& mod_divide(Number& divisor, Number& n);
		Matrix& mod_divide(Matrix& divisor, Number& n);
		
		Matrix& mod_inverse(int modulus);
		Matrix& mod_inverse(Number& modulus);
		
		Matrix& mod_pow(int exp, int m);
		Matrix& mod_pow(int exp, Number& m);
		Matrix& mod_pow(Number& exp, int m);
		Matrix& mod_pow(Number& exp, Number& m);
		
		
		
		Matrix& multiply(int val);
		Matrix& multiply(Number& number);
		Matrix& multiply(Matrix& matrix);
		Matrix& multiply(Matrix& matrix, int r);
		
		vector<Number*> multiply(vector<Number*> vec_number);
		
		
		
		Matrix& negate();
		Matrix& negate(int n);
		Matrix& negate(Number& n);
		
		Matrix& pow(int exp);
		Matrix& pow(Number& exp);
		
		int rank();
		
		int row_count();
		
		void set(Number& n, int i, int j);
		
		void set_column(vector<int> vec_int, int column);
		void set_column(vector<Number*> vec_number, int column);
		
		void set_matrix(Matrix& m, int i, int j);
		
		Matrix& set_precision(int precision);
		
		void set_row(vector<int> vec_int, int row);
		void set_row(vector<double> vec_double, int row);
		void set_row(vector<Number*> vec_number, int row);
		
		vector<Number*> solve();
		vector<Number*> solve(int n);
		vector<Number*> solve(Number& n);
		
		Matrix& square();
		
		Matrix& subtract(int n);
		Matrix& subtract(Number& n);
		Matrix& subtract(Matrix& matrix);
		
		
		Matrix& to_echelon_form();
		Matrix& to_echelon_form(int n);
		Matrix& to_echelon_form(Number& n);
		
		Matrix& to_integer();
		
		string& to_integer_string(int digits, int radix);
		
		string& to_matrix_string();
		string& to_matrix_string(int radix);
		string& to_matrix_string(int digits, int radix);
		
		Matrix& to_row_canonical_form();
		Matrix& to_row_canonical_form(Number& n);
		Matrix& to_row_canonical_form1(Number& n);
		Matrix& to_row_canonical_form2(int n);
		Matrix& to_row_canonical_form2(Number& n);
		
		string& to_string();
		string& to_string(int radix);
		string& to_string(int digits, int radix);
		
		
		Number& trace();
		
		Matrix& transpose();
		
		Matrix& trim();
	
	
	private:
	
		vector<vector<Number*>> matrix;
		
		vector<Matrix*> result;
	
	
	
	
	//  Define the overloaded assignment operator
	
	public:
	
		Matrix& operator = (Matrix &original);
	
};

//  End class Matrix

















class Math
{


	//  private Math constructor
	//
	//  no instantiation from outside of class
	
	
	private:
	
		Math();
	
	
	public:
	
	
		static constexpr double e  = 2.7182818284590452354;
		static constexpr double pi = 3.1415926535897932384;
		
		static constexpr int numberofthreads = 8;
		
		
		
		static double abs(double d);
		
		static int abs(int i);
		
		static long abs(long l);
		
		static int *add(int addend[], int elements, int summand);
		
		static vector<int> add(vector<int> a, vector<int> b);
		
		static int* add(int a[], int b[], int elements);
		
		static void add_bit(vector<int>& vec_int, long bit);
		static void add_bit(int array[], int elements, long bit);
		
		static vector<int> and_(vector<int> a1, vector<int> a2);
		static int* and_(int array1[], int array2[], int elements1, int elements2);
		
		static int bit_count(int n);
		static int bit_count(long n);
		
		static long bit_count(vector<int> vec_int);
		static long bit_count(int array[], int elements);
		
		static void clear_bit(vector<int>& vec_int, long bit);
		static void clear_bit(int array[], int elements, long bit);
		
		static int compare(vector<int> a, vector<int> b);
		static int compare(int a[], int b[], int elements1, int elements2);
		
		static int* copy(int* array, int elements);
		
		static double cos(double x);
	
	
	private:
	
		static double *cos_sin(double x);
	
	
	
	public:
	
		static double cosh(double u);
		
		static double *cos_table(int n);
		static Number& cos_table(Number& n);
		
		static double *double_array(double array[], int elements);
		
		static double e_();
		
		static bool equals(vector<int> a, vector<int> b[]);
		static bool equals(int array1[], int array2[], int elements1, int elements2);
		
		static double exp(double x);
		
		static vector<int> expand(vector<int> a, int new_size);
		
		static vector<int> factor(int n);
		static vector<int> factor(int n, int max_prime);
		static vector<int> factor(Number& n, int max_prime);
		static Number&     factor(Number& composite);
		
		static double factorial(int n);
	
	
	public:
	
		static void flip_bit(vector<int> a, long bit);
		static void flip_bit(int array[], int elements, long bit);
		
		static long gcd(long m, long n);
		
		static int get_lowest_set_bit(long a);
		
		static long get_lowest_set_bit(vector<int>& vec_int);
		static long get_lowest_set_bit(int array[], int elements);
		
		static int get_bit(vector<int>& vec_int, long bit);
		static int get_bit(int array[], int elements, long bit);
		
		static double hypot(double x, double y);
		
		static double inverse(double n);
		
		static bool is_coprime_with(long m, long n);
		
		static bool is_divisible_by(long m, long n);
		
		static bool is_power(long n, int exp);
		static bool is_power_of_2(long n);
		static bool is_power_of(int base, long n);
		
		static bool is_prime(int n);
		
		static bool is_quadratic_residue(int residue, int p);
		
		static bool is_square(double val);
		
		static bool is_square(long val);
		
		
		static int lambda(int n);
		
		static long lcm(int a, int b);
		
		static Number& lcr(vector<int> r, vector<int> n);
		static Number& lcr(int r[], int n[], int elements);
		static Number& lcr(Number r[], Number n[], int elements);
		
		
		static double log(double a, double b);
		static double log(double d);
		static double log(int a, int b);
		static double log(int d);
		static float  log(float a, float b);
		static float  log(float f);
		
		static int    log2(int n);
		static double log2(double n);
		static double log10(double n);
		
		
		static int mod(vector<int> a, int mod);
		static int mod(int array[], int elements, int mod);
		
		static int mod_divide(int dividend, int divisor, int n);
		
		static int mod_inverse(int a1, int n);
		
		static int mod_pow(int base, int exp, int mod);
		
		
		static vector<int> multiply(vector<int> a, vector<int> b);
		
		static int* multiply(int array1[], int array2[], int elements1, int elements2);
		
		static vector<int> multiply(vector<int> x, int y);
		
		static int next_prime(int n);
		
		static int next_safe_prime(int n);
		
		static Number& nCr(int n, int r);
		
		static int phi(int n);
		
		static double Poisson(double u, int x);
		
		static double pow(double base, int exp);
		static double pow(double base, long exp);
		static double pow(double base, double exp);
		
		static vector<int> primes(int n);
		
		static vector<int> quad_multiply(vector<int> x, vector<int> y);
		
		static int *quad_multiply(int x[], int y[], int elements1, int elements2);
	
	
	private:
	
		static vector<int> quad_square(vector<int> multiplicand);
		static int *quad_square(int multiplicand[], int elements);
		
		//  static Random &rng;
		
		static void init_rng();
	
	
	public:
	
		static void init_rng(long val);
		
		static double random();
		static    int random(int mod);
		static   long random(long mod);
		static double random(double mod);
		
		static double root(double n, int k);
		
		static void set_bit(vector<int>& vec_int, long bit);
		static void set_bit(int array[], int elements, long bit);
		
		static vector<int> shift_left(vector<int> a, long bits);
		static int *shift_left(int array[], int elements, long bits);
		
		static vector<int> shift_left(vector<int> a, long expansion, long bits);
		static int *shift_left(int array[], int elements, long expansion, long bits);
		
		static vector<int> shift_right(vector<int> a, long bits);
		static int *shift_right(int array[], int elements, long bits);
		
		static double *sin_table(int n);
		static Number *sin_table(Number& n);
		
		static double sin(double x);
		static double sinh(double u);
		
		
		static void sort(vector<int>& a);
		static void sort(int a[], int elements);
		static void sort(Number a[], int elements);
		static void sort(vector<vector<int>>& a);
		static void sort(int *a[], int elements);
		
		static vector<int*> sort_and_collate(int a[], int elements);
		static vector<vector<int>> sort_and_collate(vector<int> a);
		static vector<Number*> sort_and_collate(Number[], int elements);
		
		
		static double sqrt(double n);
		
		static double square(double n);
		
		static vector<int> subtract (vector<int> a, vector<int> b);
		
		static int *subtract(int minuend[], int subtrahend[], int elements);
		
		static double tan(double x);
		static double tanh(double x);
		
		static bool test_bit(vector<int> a, long bit);
		static bool test_bit(int array[], int elements, long bit);
		
		static vector<int> trim(vector<int> a);
		static int *trim(int array[], int elements);
		
		static vector<int> trim(vector<int> a, int zeros);
		static int *trim(int array[], int elements, int zeros);
		
		static vector<int> twos_complement(vector<int> a);
		static int *twos_complement(int array[], int elements);
		
		static int *unsort(int elements);
		
		static void unsort(vector<int> a);
		static void unsort(int a[], int elements);
		
		static void unsort(vector<string> *list);
};


//  End class Math










class Convert
{

	//  private Convert constructor
	//
	//  no instantiation from outside of class
	
	
	private:
	
		Convert();
	
	
	public:
	
		static const string base_16_separator;
		
		static const char int_to_base_64[64];
		static const int  base_64_to_int[128-5];
		
		
		
		static vector<char> int_vector_to_char_vector(vector<int> vec_int);
		static vector<char> char_256_vector_to_char_16_vector(vector<char> vec_char);
		
		static vector<int> int_array_to_int_vector(int array[], int elements);
		static int* int_vector_to_int_array(vector<int> vec_int);
		
		
		
		static string char_array_to_string(char array[], int elements);
		static char  *string_to_char_array(string &str);
		
		static string char_array_to_base_64(char array[], int elements);
		static char  *base_64_to_char_array(string &str);
		
		static char *char_64_array_to_char_256_array(char array[], int elements);
		static char *char_256_array_to_char_64_array(char array[], int elements);
		
		static char *char_16_array_to_char_256_array(char array[], int elements);
		static char *char_256_array_to_char_16_array(char array[], int elements);
		
		static  int *char_array_to_int_array(char array[], int elements);
		static char *int_array_to_char_array(int array[], int elements);
		
		static long *int_array_to_long_array(int array[], int elements);
		static  int *long_array_to_int_array(long array[], int elements);
		
		static double *int_array_to_complex_double_array(int array[], int elements);
		static int    *complex_double_array_to_int_array(double array[], int elements);
		
		static double *double_array_to_complex_double_array(double array[], int elements);
		static double *complex_double_array_to_double_array(double array[], int elements);
		
		static int *double_array_to_int_array(double array[], int elements);
		
		static string *string_to_string_array(string &str, int partsize);
		
		static double *double_array(double array[], int elements);
};


//  End class Convert











//  Copy control for the Matrix class
//
//  Define the copy constructor and assignment operator =


//  Copy constructor

Matrix::Matrix(Matrix& original)
{
	//  cout << "matrix copy constructor" << endl;
	
	//  vector<vector<Number*> matrix;
	//
	//  vector<Matrix*> result;
	
	for (int i = 0; i < matrix.size(); i++)
	
	    this->matrix[i] = original.matrix[i];
}



//  Assignment operator

Matrix& Matrix::operator = (Matrix& original)
{
	//  cout << "matrix assignment operator = " << endl;
	
	//  vector<vector<Number*> matrix;
	//
	//  vector<Matrix*> result;
	
	for (int i = 0; i < matrix.size(); i++)
	
	    this->matrix[i] = original.matrix[i];
	
	return *this;
}







//  Matrix destructor
//
//
//  The Matrix destructor is called recursively for equations such as
//
//  A = A .function() .function() .function() .function() ...
//
//  because it creates pointers to pointers to pointers ... to Matrix,
//  or a linked list of pointers to Matrix.
//
//  As each function is evaluated, it returns a new anonymous object
//  which is used to evaluate the next function, which returns another
//  new anonymous object which is then used to evaluate the next function,
//  ... until the last new object is assigned to the variable.



Matrix::~Matrix()
{
	//  vector<vector<Number*> matrix;
	//
	//  vector<Matrix*> result;
	
	//  cout << "entering matrix destructor" << endl;
	
	//  Delete the matrix of Number pointers
	
	for (int i = 0; i < this->matrix   .size(); i++)
	for (int j = 0; j < this->matrix[i].size(); j++)
	
	    delete this->matrix[i][j];
	
	//  Delete the pointers to pointers to matrix
	//     and then delete the pointers to matrix
	
	for (int i = 0; i < this->result.size(); i++)
	{
		if (this->result[i] != nullptr)
		{
			for (int j = 0; j < this->result[i]->result.size(); j++)
			{
				if (this->result[i]->result[j] != nullptr)
				{
					//  cout << "deleting pointer to pointer" << endl;
					
					delete this->result[i]->result[j];
			        	       this->result[i]->result[j] = nullptr;
				}
			}
			
			delete this->result[i];
			       this->result[i] = nullptr;
		}
	}
	
	//  cout << "exiting matrix destructor" << endl << endl;
}




//  Matrix constructors



Matrix::Matrix()
{

}


Matrix::Matrix(int rows, int columns)
{

	this->matrix = vector<vector<Number*>> (rows);
	
	//  Initialize each row to the number of columns
	
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			Number *zero = new Number(0);
			
			this->matrix[i].push_back(zero);
		}
	}
}



Matrix::Matrix(vector<int> matrix, int rows, int columns)
{

	this->matrix = vector<vector<Number*>> (rows);
	
	//  Initialize each row to the number of columns
	
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			int index = columns * i + j;
			
			Number *element = new Number(matrix[index]);
			
			this->matrix[i].push_back(element);
		}
	}
}


Matrix::Matrix(vector<double> matrix, int rows, int columns)
{

	//  ...
	
	//  ...
}


Matrix::Matrix(vector<Number*> matrix, int rows, int columns)
{

	//  ...
	
	//  ...
}



//  This constructor creates a horizontal matrix / single row
//
//  Use  Matrix(array) .transpose() to get a vertical matrix / single column


Matrix::Matrix(vector<Number*> row)
{

	//  ...
	
	//  ...
}


Matrix::Matrix(vector<int> row)
{

	//  ...
	
	//  ...
}


Matrix::Matrix(vector<double> row)
{

	//  ...
	
	//  ...
}



Matrix::Matrix(string &str, int radix, int rows, int columns)
{

	//  This constructor is used in cryptography to convert
	//
	//  a public key string to a public key matrix
	
	
	//  Initialize the matrix to the number of rows and
	//  initialize the rows to the number of columns
	
	this->matrix = vector<vector<Number*>>(rows);
	
	for (int i = 0; i < rows; i++)
	
	    this->matrix[i] = vector<Number*>(columns);
	
	
	int pos = 0;
	
	while ((pos = str.find_first_of(" ")) != string::npos)
	
	    str = str .erase(pos, 1);
	
	if (   ( str.length() < (rows * columns) )
	  || ( ( str.length() % (rows * columns) ) != 0 ) )
	{
		string message = "matrix string length"
		
		    + string(" is not a multiple of rows x columns");
	
		cout << message << endl; throw string(message);
	}
	
	int digits = str.length() / (rows * columns);
	
	//  For each row copy a substring of string length / rows
	
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			string substring = str.substr(
			
			    (i * columns + j + 0) * digits,
			    (i * columns + j + 1) * digits);
			
			Number number(substring, radix);
			
			this->matrix[i][j] = new Number(number);
		}
	}
}






//  Matrix methods



Matrix& Matrix::set_result(Matrix& matrix)
{

	Matrix *result = new Matrix(matrix);
	
	this->result.push_back(result);
	
	return *result;
}





Matrix& Matrix::abs()
{
	//  returns the absolute value
	
	Matrix matrix (*this);
	
	for (int i = 0; i < matrix.matrix   .size(); i++)
	for (int j = 0; j < matrix.matrix[i].size(); j++)
	
	    matrix.matrix[i][j] = &
	    matrix.matrix[i][j]->abs();
	
	return set_result(matrix);
}


Matrix& Matrix::add(int n)
{
	Number number(n);
	
	return this->add(n);
}


Matrix& Matrix::add(Number& n)
{
	//  adds a number n to a matrix
	
	Matrix matrix(*this);
	
	for (int i = 0; i < matrix.matrix   .size(); i++)
	for (int j = 0; j < matrix.matrix[i].size(); j++)
	
	    matrix.matrix[i][j] = &
	    matrix.matrix[i][j] ->add(n);
	
	return set_result(matrix);
}


Matrix& Matrix::add(Matrix& matrix)
{
	//  adds two matrices
	
	Matrix matrix1(*this);
	
	for (int i = 0; i < matrix1.matrix   .size(); i++)
	for (int j = 0; j < matrix1.matrix[i].size(); j++)
	
	    matrix1.matrix[i][j] = &
	    matrix1.matrix[i][j]
	
		-> add(*matrix.matrix[i][j]);
	
	return set_result(matrix1);
}


Matrix& Matrix::augment(vector<int> b)
{
	//  augments a matrix
	
	int size = b.size();
	
	vector<Number*> B (size);
	
	for (int i = 0; i < size; i++)
	
	    B[i] = new Number(b[i]);
	
	return this->augment(B);
}


Matrix& Matrix::augment(vector<Number*> B)
{

	//  augments a matrix
	
	Matrix A (*this);
	
	if (B.size() != A.matrix.size())
	{
		string message =
		
		  "length of number != number of rows in matrix";
		
		cout << message << endl; throw string(message);
	}
	
	//  Expand the matrix by one column
	
	Matrix matrix (A.matrix.size(), A.matrix[0].size() + 1);
	
	for (int i = 0; i < A.matrix   .size(); i++)
	for (int j = 0; j < A.matrix[i].size(); j++)
	
	    matrix.matrix[i][j] = A.matrix[i][j];
	
	
	//  Set the last column equal to B
	
	for (int i = 0; i < matrix.matrix.size(); i++)
	
	    matrix.matrix[i][matrix.matrix[i].size() -1] = B[i];
	
	return set_result(matrix);
}


Matrix& Matrix::augment(Matrix& B)
{

	//  augments a matrix
	
	Matrix A (*this);
	
	if (A.matrix.size() != B.matrix.size())
	{
		string message = "number of rows in matrix B != A";
		
		cout << message << endl; throw string(message);
	}
	
	Matrix matrix (A.row_count(), A.column_count() + B.column_count());
	
	for (int i = 0; i < A.matrix   .size(); i++)
	for (int j = 0; j < A.matrix[i].size(); j++)
	
	    matrix.matrix[i][j] = A.matrix[i][j];
	
	for (int i = 0; i < B.matrix   .size(); i++)
	for (int j = 0; j < B.matrix[i].size(); j++)
	
	    matrix.matrix[i][A.matrix[i].size() + j] = B.matrix[i][j];
	
	return set_result(matrix);
}


int Matrix::column_count()
{
	if (this->matrix.size() == 0) return 0;
	
	return this->matrix[0].size();
}


Matrix& Matrix::delete_column(int column)
{

	int rows = this->row_count();
	
	int columns = this->column_count();
	
	Matrix matrix (rows, columns -1);
	
	for (int i = 0;        i < rows;       i++)
	for (int j = 0, k = 0; j < columns -1; j++, k++)
	{
		if (j == column)  k++;
		
		matrix.matrix[i][j] = this->matrix[i][k];
	}
	
	return set_result(matrix);
}


Matrix& Matrix::delete_row(int row)
{

	//  Delete row is used by the determinant method
	
	//  and by the trim method
	
	int rows = this->row_count();
	
	int columns = this->column_count();
	
	Matrix matrix (rows -1, columns);
	
	for (int this_i = 0, matrix_i = 0; this_i < rows; this_i++)
	
	    if (this_i != row) { matrix.set_row(
	
		this->get_row(this_i), matrix_i); matrix_i++; }
	
	return set_result(matrix);
}


Number& Matrix::determinant()
{
	Number zero(0);
	
	return this->determinant(zero);
}


Number& Matrix::determinant(Number& n)
{

	//  computes the determinant of a matrix by reducing the
	//  matrix to echelon (upper triangular) form and then
	//  computing the product of the diagonal elements
	
	
	if (this->row_count() != this->column_count())
	{
		string message = "non-square matrix determinant";
		
		cout << message << endl; throw string(message);
	}
	
	int swaps = 0;
	
	Number determinant(1);
	
	
	//  Reduce the matrix to echelon form
	
	Matrix matrix(*this);
	
	int precision = Number(1) .inverse() .get_precision();
	
	for (int j = 0; j < matrix.column_count(); j++)
	
	    matrix.matrix[0][j] = & matrix.matrix[0][j]
	
		->set_precision(precision);
	
	for (int r = 0; r < matrix.matrix.size() -1; r++)
	{
		//  Find the row and column of the first non-zero
		//  element starting from row i = r, column j = r
		
		bool bool1 = false;
		
		int i = 0, j = 0;
		
		for (j = r; j < matrix.column_count(); j++)
		{
			for (i = r; i < matrix.row_count(); i++)
			{
				if (!n.equals(0))
				{
					if (!matrix.matrix[i][j]->equals(0))
					
					    { bool1 = true;  break; }
				}
				
				else if (n.equals(0))
				{
					if (!matrix.matrix[i][j]->equals(0.0))
					
					    { bool1 = true;  break; }
				}
			}
			
			if (bool1 == true) break;
		}
		
		if (j >= matrix.column_count()) continue;
		
		
		//  Swap rows i and r if element [r][j] equals zero
		
		if (i != r)
		{
			swaps += 1;
			
			vector<Number*> tempi = matrix.matrix[i];
			vector<Number*> tempr = matrix.matrix[r];
			
			matrix.matrix[i] = tempr;
			matrix.matrix[r] = tempi;
		}
		
		
		//  Put zeros below the pivot for each row below r
		
		for (i = r + 1; i < matrix.row_count(); i++)
		{
			//  Compute the multiplier m for row i
			
			Number m, m1, m2;
			
			if (n.equals(0))
			{
				m = matrix.matrix[i][j] ->negate() .multiply(
				
				    matrix.matrix[r][j] ->inverse() );
			}
			
			else if (!n.equals(0))
			{
				m1 =   matrix.matrix[i][j] ->negate(n);
				
				m2 = * matrix.matrix[r][j];
			}
			
			
			//  Multiply row r by the multiplier, and add to row i
			
			for (int k = j; k < matrix.column_count(); k++)
			{
				if (n.equals(0))
				{
					matrix.matrix[i][k] = & matrix.matrix[i][k]
					
					    ->add( m.multiply( *matrix.matrix[r][k] ) );
				}
				
				else if (!n.equals(0))
				{
					matrix.matrix[i][k] = & m1.multiply(*matrix.matrix[r][k])
					
					    .add( m2.multiply( *matrix.matrix[i][k] ) )
					
						.mod(n).add(n).mod(n);
					
					determinant = determinant .multiply(m2) .mod(n).add(n).mod(n);
				}
			}
		}
	}
	
	
	//  Calculate the product of the diagonal elements
	
	for (int k = 0; k < matrix.row_count(); k++)
	{
		determinant = determinant .multiply(*matrix.matrix[k][k]);
		
		if (!n.equals(0)) determinant = determinant.mod(n);
	}
	
	if ((swaps % 2) == 1)
	
	    determinant = determinant.negate();
	
	return *new Number(determinant);
}


Number& Matrix::determinant1()
{
	Number zero(0);
	
	return this->determinant1(zero);
}

Number& Matrix::determinant1(Number& n)
{

	//  computes the determinant of a matrix
	//  by expanding the matrix into minors
	
	//  Compute the determinant by recursion
	//
	//  using the definitions | a | == a
	//
	//      | a  b |
	//  and |      | == a d - b c.
	//      | c  d |
	
	
	if (!this->is_square())
	{
		string message = "non-square matrix";
		
		cout << message << endl; throw string(message);
	}
	
	Matrix matrix(*this);
	
	if (!n.equals(0)) matrix = matrix.mod(n);
	
	
	//  Test before recursion
	
	if (matrix.row_count() == 1)
	
	    return *matrix.matrix[0][0];
	
	
	if (matrix.row_count() == 2)
	{
		//  det == a[0][0] a[1][1] - a[0][1] a[1][0]
		
		Number determinant = matrix.matrix[0][0]
		
		    ->multiply(*matrix.matrix[1][1])
		
		        .subtract(matrix.matrix[0][1]
		
			    ->multiply(*matrix.matrix[1][0]));
		
		if (!n.equals(0))
		
		    determinant = determinant .mod(n).add(n).mod(n);
		
		return *new Number(determinant);
	}
	
	
	//  Start recursion
	
	//  determinant = the sum of A[0][i] * the
	//
	//     determinant of A.delete_row(0) .delete_column(i)
	
	
	Number determinant(0);
	
	for (int i = 0; i < this->column_count(); i++)
	{
		Number determinant1 = matrix.matrix[0][i] ->multiply(
		
		    matrix .delete_row(0) .delete_column(i) .determinant1(n) );
		
		if ((i % 2) == 1)  determinant1 = determinant1 .negate();
		
		determinant = determinant .add(determinant1);
	}
	
	return *new Number(determinant);
}


Matrix& Matrix::divide(int n)
{
	Number number(n);
	
	return this->divide(number);
}


Matrix& Matrix::divide(Number& number)
{

	Matrix matrix(*this);
	
	for (int i = 0; i < matrix.row_count(); i++)
	for (int j = 0; j < matrix.column_count(); j++)
	
	    matrix.matrix[i][j] = & this->matrix[i][j] ->divide(number);
	
	return set_result(matrix);
}


bool Matrix::equals(Matrix& matrix)
{
	if ( (this->   row_count() != matrix.   row_count())
	  || (this->column_count() != matrix.column_count()) )
	{
		string message = "matrices are not comparable";
		
		cout << message << endl; throw string(message);
	}
	
	for (int i = 0; i < this->   row_count(); i++)
	for (int j = 0; j < this->column_count(); j++)
	
	    if (! this->matrix[i][j]->equals(*matrix.matrix[i][j]))
	
		return false;
	
	return true;
}


Number& Matrix::get(int i, int j)
{
	//  returns one element of a matrix
	
	return *this->matrix[i][j];
}


Matrix& Matrix::get(int row, int column, int rows, int columns)
{
	//  copies a sub matrix from a (super) matrix
	
	//  Row and column are the offsets from [0, 0]
	
	if ( (   row + rows    > this->   row_count())
	  || (column + columns > this->column_count()) )
	{
		string message = "rows or columns of sub matrix + offset"
		        + string(" exceed the super matrix dimensions");
		
		cout << message << endl; throw string(message);
	}
	
	Matrix matrix (rows, columns);
	
	for (int i = 0; i < matrix.   row_count(); i++)
	for (int j = 0; j < matrix.column_count(); j++)
	
	    matrix.set(*this->matrix[row + i][column + j], i, j);
	
	return set_result(matrix);
}


vector<Number*> Matrix::get_column(int column)
{
	//  returns one column of a matrix
	
	vector<Number*> vec_number (row_count());
	
	for (int i = 0; i < row_count(); i++)
	
	    vec_number[i] = matrix[i][column];
	
	return vec_number;
}


vector<Number*> Matrix::get_diagonals()
{
	//  returns the diagonal components
	
	vector<Number*> diagonals (row_count());
	
	for (int i = 0; i < row_count(); i++)
	
	    diagonals[i] = matrix[i][i];
	
	return diagonals;
}


int Matrix::get_precision()
{
	returns the precision
	
	int precision = 0;
	
	for (int i = 0; i < row_count(); i++)
	for (int j = 0; j < column_count(); j++)
	
	    if (matrix[i][j]->get_precision() > precision)
	
		precision = matrix[i][j]->get_precision();
	
	return precision;
}


vector<Number*> Matrix::get_row(int row)
{
	//  returns one row of a matrix
	
	int columns = this->column_count();
	
	vector<Number*> vec_number (columns);
	
	for (int j = 0; j < columns; j++)
	
	    vec_number[j] = matrix[row][j];
	
	return vec_number;
}



Matrix& Matrix::identity_matrix(int size)
{
	//  returns the unit matrix
	//
	//  | 1                          |
	//  |    1                       |
	//  |       1                    |
	//  |          1                 |
	//  |             1              |
	//  |                1           |
	//  |                   .        |
	//  |                      .     |
	//  |                         .  |
	
	Matrix matrix (size, size);
	
	for (int i = 0; i < size; i++)
	for (int j = 0; j < size; j++)
	{
		if (j != i)
		
		     matrix.matrix[i][j] = new Number(0);
		else matrix.matrix[i][j] = new Number(1);
	}
	
	return set_result(matrix);
}


Matrix& Matrix::inverse()
{

	//  computes the inverse of a matrix
	
	//  To solve the equation A X == I, form the augmented matrix
	//  M = [ A , I ] and solve the system by reducing the matrix
	//  M to echelon form and then to row canonical form.
	//
	//  The left side will be the identity matrix I and
	//  the right side will be the inverse matrix A ^-1.
	
	
	Matrix matrix(*this);
	
	int rows    = matrix.row_count();
	int columns = matrix.column_count();
	
	if (!matrix .is_square())
	{
		string message = "non-square matrix inversion";
		
		cout << message << endl; throw string(message);
	}
	
	
	//  Augment the matrix with the identity matrix
	
	Matrix I = identity_matrix(rows);
	
	matrix = matrix .augment(I);
	
	
	//  Use Gaussian elimination to invert the matrix
	
	
	//  Forward elimination
	//
	//  Reduce the matrix to upper echelon form
	//
	//  (put zeros below each pivot element)
	
	matrix = matrix.to_echelon_form();
	
	
	//  Back substitution
	//
	//  (put zeros above each pivot element)
	
	matrix = matrix.to_row_canonical_form();
	
	
	//  The right half of the matrix contains the inverse of M
	
	Matrix inv_matrix = matrix.get(0, columns, rows, columns);
	
	return set_result(inv_matrix);
}


bool Matrix::is_complex()
{
	//  tests if the matrix is complex
	
	for (int i = 0; i < this->row_count(); i++)
	for (int j = 0; j < this->column_count(); j++)
	
	    if (matrix[i][j]->is_complex()) return true;
	
	return false;
}


bool Matrix::is_echelon_form()
{
	//  tests if the matrix is in echelon form
	
	//  If a matrix is in echelon form, the first non-zero element
	//  of each row is to the right of the first non-zero element
	//  of the preceding row.
	
	Matrix matrix(*this);
	
	Number zero(0);
	
	if (matrix.get_precision() != 0)
	
	    zero = zero .set_precision(matrix.get_precision());
	
	int index = -1;
	
	for (int i = 0; i < matrix.row_count();       i++)
	for (int j = 0; j < matrix.column_count() -1; j++)
	{
		if (!matrix.matrix[i][j] ->equals(zero))
		{
			if (j <= index) return false;
			
			index = j;  break;
		}
	}
	
	return true;
}


bool Matrix::is_identity_matrix()
{

	Number  one(1);
	Number zero(0);
	
	bool is_integer = matrix[0][0]->is_integer();
	
	if (!is_integer)
	{
		int precision = matrix[0][0]->get_precision();
		
		one  = one  .set_precision(precision);
		zero = zero .set_precision(precision);
	}
	
	for (int i = 0; i < this->row_count(); i++)
	for (int j = 0; j < this->column_count(); j++)
	{
	    if      ((i == j) && !matrix[i][j]->equals(one))  return false;
	    else if ((i != j) && !matrix[i][j]->equals(zero)) return false;
	}
	
	return true;
}



bool Matrix::is_row_canonical_form()
{

	//  tests if the matrix is in row canonical form
	
	//  If a matrix is in row canonical form, the diagonal elements
	//  have to equal 1 and the non-diagonal elements have to equal
	//  0 except for the last column or right side of the matrix.
	
	Matrix matrix(*this);
	
	Number zero (0);
	Number  one (1);
	
	if (matrix.get_precision() != 0)
	{
	    zero = zero .set_precision(matrix.get_precision());
	    one  = one  .set_precision(matrix.get_precision());
	}
	
	for (int i = 0; i < matrix.column_count() -1; i++)
	
	    if (!matrix.matrix[i][i]->equals(one)) return false;
	
	for (int i = 0; i < matrix.column_count() -1; i++)
	for (int j = 0; j < matrix.column_count() -1; j++)
	
	    if ((i != j) && !matrix.matrix[i][j] ->equals(zero))
	
		return false;
	
	return true;
}


bool Matrix::is_singular()
{
	int p = this->get_precision();
	
	Number zero(0); zero = zero .set_precision(p);
	
	return this->determinant().abs().equals(zero);
}

bool Matrix::is_singular(Number& n)
{
	return this->determinant(n).equals(0);
}

bool Matrix::is_square()
{
	return row_count() == column_count() ? true : false;
}

Matrix *Matrix::LU()
{
	Number zero(0);
	
	return this->LU(zero);
}


Matrix *Matrix::LU(Number& modulus)
{

	//  decomposes a square matrix into lower and
	//  upper triangular matrices
	
	//  The upper triangular matrix U is computed from
	//  a matrix A by reducing A to echelon form.
	//
	//  The lower triangular matrix L is defined
	//  (on, above, and below the diagonal) by
	//
	//  L[i][j] =   1          i == j,
	//  L[i][j] =   0          i <  j,
	//  L[i][j] = - m[i][j]    i >  j,
	//
	//  where m[i][j] are the multipliers used to reduce the
	//  matrix A to echelon or upper triangular form U.
	
	
	Number n = modulus;
	
	Matrix A(*this);
	
	if (A.row_count() != A.column_count())
	{
		string message = "LU factorization of non-square matrix";
		
		cout << message << endl; throw string(message);
	}
	
	int s = A.row_count();
	
	Number m[s][s];
	
	Matrix L(s, s), U(*this);
	
	
	//  Compute the upper matrix
	
	for (int r = 0; r < U.row_count(); r++)
	{
		//  Find the row and column of the first non-zero element
		//  starting from row i = r and column j = r.
		
		bool bool1 = false;
		
		int i = 0, j = 0;
		
		for (j = r; j < U.column_count(); j++)
		{
			//  columns
			
			for (i = r; i < U.row_count(); i++)
			{
				//  rows
				
				if (!U.matrix[i][j]->equals(0))
				
				    { bool1 = true; break; }
			}
			
			if (bool1 == true) break;
		}
		
		if (j == U.column_count()) continue;
		
		
		//  Swap rows i and r if element [r][j] is zero
		
		if (i != r)
		{
			vector<Number*> tempi = U.matrix[i];
			vector<Number*> tempr = U.matrix[r];
			
			U.matrix[i] = tempr;
			U.matrix[r] = tempi;
		}
		
		
		//  Put zeros below the pivot for each row below r
		
		for (i = r+1; i < U.row_count(); i++)
		{
			//  Compute the multiplier for row i
			
			m[i][j] = U .matrix[i][j] ->negate()
			
			    .multiply(((!n.equals(0)) ?
			
				U.matrix[r][j] ->mod_inverse(n) :
				U.matrix[r][j] ->    inverse()));
			
			
			//  Multiply row i by the multiplier, and add to row i
			
			for (int k = j; k < (U.column_count()); k++)
			{
				U.matrix[i][k] = & U.matrix[i][k]
				
				    ->add( m[i][j] .multiply(*U.matrix[r][k]) );
				
				if (!n.equals(0))
				
				    U.matrix[i][k] = & U.matrix[i][k]
				
					->mod(n).add(n).mod(n);
			}
		}
	}
	
	
	//  Compute the lower matrix
	
	for (int i = 0; i < s; i++)
	for (int j = 0; j < s; j++)
	{
		if (i > j)
		{
			L.matrix[i][j] = & m[i][j].negate();
			
			if (!n.equals(0)) L.matrix[i][j] = & L.matrix[i][j]
			
			    ->mod(n).add(n).mod(n);
		}
		
		else if (i == j) L.matrix[i][j] = new Number(1);
		else if (i <  j) L.matrix[i][j] = new Number(0);
	}
	
	Matrix *LU = new Matrix[2];
	
	LU[0] = L; LU[1] = U;
	
	return LU;
}


Matrix& Matrix::mod(int n)
{
	Number number(n);
	
	return this->mod(number);
}

Matrix& Matrix::mod(Number& modulus)
{
	Matrix matrix(*this);
	
	for (int i = 0; i < matrix.row_count(); i++)
	for (int j = 0; j < matrix.column_count(); j++)
	
	    matrix.matrix[i][j] = & this->matrix[i][j] ->mod(modulus);
	
	return set_result(matrix);
}


Matrix& Matrix::mod_divide(int divisor, Number& mod)
{
	Number number(divisor);
	
	return mod_divide(number, mod);
}


Matrix& Matrix::mod_divide(Number& divisor, Number& n)
{
	Matrix a = this-> mod(n).add(n).mod(n);
	
	if (!divisor.is_integer() || !n.is_integer())
	{
		string message = "non-integer divisor or modulus";
		
		cout << message << endl; throw string(message);
	}
	
	if (!divisor.is_coprime_with(n))
	{
		string message = "divisor is non-invertible";
		
		cout << message << endl; throw string(message);
	}
	
	Matrix quotient = a .multiply( divisor.mod_inverse(n) ) .mod(n);
	
	if (! quotient .multiply(divisor) .mod(n).add(n).mod(n) .equals(a) )
	{
		string message = "matrix mod division error";
		
		cout << message << endl; throw string(message);
	}
	
	return set_result(quotient);
}


Matrix& Matrix::mod_divide(Matrix& divisor, Number& n)
{
	Matrix a = this-> mod(n).add(n).mod(n);
	
	Matrix quotient = a .multiply(divisor.mod_inverse(n)) .mod(n);
	
	if (! quotient .multiply(divisor) .mod(n).add(n).mod(n) .equals(a) )
	{
		string message = "matrix modular division error";
		
		cout << message << endl; throw string(message);
	}
	
	return set_result(quotient);
}


Matrix& Matrix::mod_inverse(int modulus)
{
	Number number(modulus);
	
	return mod_inverse(number);
}

Matrix& Matrix::mod_inverse(Number& modulus)
{

	//  computes the inverse of a matrix A
	//  or solves the equation A X == I (mod n)
	
	Number n = modulus;
	
	Matrix matrix(*this);
	
	matrix = matrix.mod(n);
	
	if (!matrix .is_square())
	{
		string message = "non-square matrix inversion";
		
		cout << message << endl; throw string(message);
	}
	
	int rows    = matrix.row_count();
	int columns = matrix.column_count();
	
	
	if (this->determinant(n).equals(0))
	{
		string message = "non-invertible matrix";
		
		message += " determinant == 0";
		
		cout << message << endl; throw string(message);
	}
	
	
	//  Augment the matrix with the identity matrix
	
	Matrix I = identity_matrix(rows);
	
	matrix = matrix .augment(I);
	
	
	//  Use Gaussian elimination to invert the matrix
	
	//  Forward elimination (upper triangularization)
	
	//  Reduce the matrix to (upper) echelon form
	//  (put zeros below each pivot variable)
	
	if (modulus.trim().length() > 1)
	
	     matrix = matrix.to_echelon_form(n);
	else matrix = matrix.to_echelon_form(n.int_value());
	
	
	//  Back substitution (lower triangularization)
	
	//  Reduce the matrix to row canonical form
	//  (put zeros above each pivot variable)
	
	matrix = matrix.to_row_canonical_form(n);
	
	
	//  Copy the inverse of M from the right half of the augmented matrix
	//
	//  (The inverse matrix is the same size as the original matrix)
	
	Matrix inv_matrix = matrix.get(0, columns, rows, columns);
	
	
	//  Verify the modular inverse
	
	if (!this->multiply(inv_matrix) .mod(n) .add(n) .mod(n) .equals(I))
	{
		string message = "matrix mod inverse error";
		
		cout << message << endl; throw string(message);
	}
	
	return set_result(inv_matrix);
}



Matrix& Matrix::mod_pow(int exp, int m)
{
	Number n1 (exp), n2 (m);
	
	return mod_pow(n1, n2);
}

Matrix& Matrix::mod_pow(int exp, Number& m)
{
	Number n1 (exp), n2 (m);
	
	return mod_pow(n1, n2);
}

Matrix& Matrix::mod_pow(Number& exp, int m)
{
	Number n1 (exp), n2 (m);
	
	return mod_pow(n1, n2);
}

Matrix& Matrix::mod_pow(Number& exp, Number& m)
{
	//  computes A ^ x (mod m)
	
	if (this->row_count() != this->column_count())
	{
		string message = "non-square matrix exponentiation";
		
		cout << message << endl; throw string(message);
	}
	
	Matrix a (*this); a = a.mod(m);
	
	Number x (exp);
	
	Matrix y = identity_matrix(this->row_count());
	
	
	//  Compute y = a ^ x  (mod m)
	
	//  Use the square and multiply method
	
	while (!x.equals(0))
	{
		//  Accumulate squares
		
		if (x.test_bit(0))
		
		    y = a.multiply(y) .mod(m);
		
		//  Square the square
		
		a = a.square(). mod(m);
		
		//  Shift the exponent to the next bit
		
		x = x.shift_right(1);
	}
	
	return set_result(y);
}


Matrix& Matrix::multiply(int val)
{
	//  multiplies a matrix by an int
	
	Number number(val);
	
	return this->multiply(number);
}

Matrix& Matrix::multiply(Number& number)
{
	//  multiplies a matrix by a number
	
	Matrix matrix(*this);
	
	for (int i = 0; i < matrix.row_count(); i++)
	for (int j = 0; j < matrix.column_count(); j++)
	
	    matrix.matrix[i][j] = &
	    matrix.matrix[i][j] ->multiply(number);
	
	return set_result(matrix);
}

Matrix& Matrix::multiply(Matrix& matrix)
{
	//  multiplies a matrix by a matrix
	
	if (this->column_count() != matrix.row_count())
	{
		string message = "number of columns of first matrix "
		
		    + string("does not equal rows of second matrix");
		
		cout << message << endl; throw string(message);
	}
	
	//  Set the dimensions of the product matrix
	
	Matrix product (this->row_count(), matrix.column_count());
	
	//  Compute the product elements
	//
	//  c[i, j] = the sum of  a[i, k] b[k, j]
	
	for (int i = 0; i < product   .row_count(); i++)
	for (int j = 0; j < product.column_count(); j++)
	for (int k = 0; k <   this->column_count(); k++)
	
	    product.matrix[i][j] = & product.matrix[i][j]
	
		->add(this->matrix[i][k] ->multiply(
		    *matrix.matrix[k][j]));
	
	return set_result(product);
}


vector<Number*> Matrix::multiply(vector<Number*> X)
{
	//  multiplies a matrix by a vector
	//
	//  This is used to verify the equation A X == B.
	
	if (this->column_count() != X.size())
	{
		string message = "number of columns of first matrix "
		
		    + string("does not equal rows of second matrix");
		
		cout << message << endl; throw string(message);
	}
	
	//  Multiply each row of A by column X and add the products to calculate B[i]
	
	//  The length of the product (column) equals the number of matrix rows
	
	int rows = this->row_count();
	
	int columns = this->column_count();
	
	vector<Number*> B (rows);
	
	vector<Number*> temp (rows);
	
	for (int i = 0; i < rows; i++)
	
	    temp[i] = new Number(0);
	
	for (int i = 0; i < rows; i++)
	
	   for (int j = 0; j < columns; j++)
	
	      temp[i] = & temp[i] ->add(
	
		 this->matrix[i][j]->multiply(*X[j]) );
	
	for (int i = 0; i < rows; i++)
	
	    B[i] = temp[i];
	
	return B;
}



Matrix& Matrix::multiply(Matrix& matrix, int r)
{
	//  multiplies the first r rows of a matrix by a matrix
	
	if (this->column_count() != matrix.row_count())
	{
		string message = "number of columns of first matrix "
		
		  + string("does not equal rows of second matrix");
		
		cout << message << endl; throw string(message);
	}
	
	Matrix product (row_count(), column_count());
	
	for (int i = 0; i < r; i++)
	
	for (int j = 0; j < product.column_count(); j++)
	for (int k = 0; k <   this->column_count(); k++)
	
	    product.matrix[i][j] = & product.matrix[i][j] ->add(
	
	      this->matrix[i][k]->multiply(*matrix.matrix[k][j]));
	
	return set_result(product);
}



Matrix& Matrix::negate()
{
	Matrix matrix(*this);
	
	for (int i = 0; i < matrix.row_count(); i++)
	for (int j = 0; j < matrix.column_count(); j++)
	
	    matrix.matrix[i][j] = &
	    matrix.matrix[i][j]->negate();
	
	return set_result(matrix);
}

Matrix& Matrix::negate(int n)
{
	Number number(n);
	
	return this->negate(number);
}


Matrix& Matrix::negate(Number& n)
{
	Matrix matrix(*this);
	
	for (int i = 0; i < matrix.row_count(); i++)
	for (int j = 0; j < matrix.column_count(); j++)
	
	    matrix .matrix[i][j] = &
	    matrix .matrix[i][j]->negate(n);
	
	return set_result(matrix);
}


Matrix& Matrix::pow(int exp)
{
	Number number(exp);
	
	return this->pow(number);
}


Matrix& Matrix::pow(Number& exp)
{
	if (this->row_count() != this->column_count())
	{
		string message = "non-square matrix exponentiation";
		
		cout << message << endl; throw string(message);
	}
	
	//  Define the square Matrix a and exponent x
	
	Matrix a (*this);
	
	Number x (exp);
	
	
	//  Define the output matrix y and initialize to one
	
	Matrix y = identity_matrix(this->row_count());
	
	
	//  Compute y = a^x
	
	//  Use the square and multiply method
	
	while (!x.equals(0))
	{
		//  Accumulate squares
		
		if (x.test_bit(0))
		
		    y = a.multiply(y);
		
		//  Square the square
		
		a = a.square();
		
		//  Shift the exponent to the next bit
		
		x = x.shift_right(1);
	}
	
	return set_result(y);
}


int Matrix::rank()
{
	//  returns the rank of a matrix
	
	//  The rank of a matrix is the number of linearly independent equations
	
	return to_echelon_form() .trim() .row_count();
}


int Matrix::row_count()
{
	return this->row_count();
}


void Matrix::set(Number& n, int i, int j)
{
	if (this->matrix[i][j] != nullptr)
	
	    delete this->matrix[i][j];
	
	this->matrix[i][j] = new Number(n);
}

void Matrix::set_column(vector<int> vec_int, int column)
{
	int size = vec_int.size();
	
	if (size != this->row_count())
	{
		string message = "vector length != number of rows";
		
		cout << message << endl; throw string(message);
	}
	
	vector<Number*> vec_number (size);
	
	for (int i = 0; i < size; i++)
	
	    vec_number[i] = new Number(vec_int[i]);
	
	this->set_column(vec_number, column);
}


void Matrix::set_column(vector<Number*> vec_number, int column)
{
	int j = column;
	
	for (int i = 0; i < this->row_count(); i++)
	
	    this->matrix[i][j] = vec_number[i];
}

void Matrix::set_matrix(Matrix& m, int i, int j)
{
	//  pastes a smaller matrix into a larger
	//  matrix starting at (i, j)
	
	//  The size of the sub matrix plus the offset
	//  cannot exceed the size of the super matrix
	
	if (((i + m.row_count())    > this->row_count())
	 || ((j + m.column_count()) > this->column_count()))
	{
		string message = "rows or columns of sub matrix + offset"
		
		    + string(" exceed dimensions of super matrix");
		
		cout << message << endl; throw string(message);
	}
	
	for (int r = 0; r < m.row_count(); r++)
	for (int s = 0; s < m.column_count(); s++)
	
	    this->matrix[i+r][j+s] = m.matrix[r][s];
}


Matrix& Matrix::set_precision(int precision)
{
	Matrix matrix(*this);
	
	for (int i = 0; i < matrix.row_count(); i++)
	for (int j = 0; j < matrix.column_count(); j++)
	
	    matrix.matrix[i][j] = &
	    matrix.matrix[i][j]
	
		->set_precision(precision);
	
	return set_result(matrix);
}


void Matrix::set_row(vector<int> vec_int, int row)
{
	int size = vec_int.size();
	
	vector<Number*> vec_number (size);
	
	for (int i = 0; i < size; i++)
	
	    vec_number[i] = new Number(vec_int[i]);
	
	this->set_row(vec_number, row);
}

void Matrix::set_row(vector<double> vec_double, int row)
{
	int size = vec_double.size();
	
	vector<Number*> vec_number (size);
	
	for (int i = 0; i < size; i++)
	
	    vec_number[i] = new Number(vec_double[i]);
	
	this->set_row(vec_number, row);
}

void Matrix::set_row(vector<Number*> vec_number, int row)
{
	int size = vec_number.size();
	
	vector<Number*> vector1 (size);
	
	for (int i = 0; i < size; i++)
	
	    vector1[i] = vec_number[i];
	
	if (!this->matrix[row].empty())
	{
		//  Make sure the caller is not making a mistake by resetting a row
		
		string message = "initialize row to nullptr before setting row";
		
		cout << message << endl; throw string(message);
	}
	
	this->matrix[row] = vector1;
}



vector<Number*> Matrix::solve()
{
	Number mod(0);
	
	return this->solve(mod);
}


vector<Number*> Matrix::solve(int n)
{
	Number number(n);
	
	return this->solve(number);
}


vector<Number*> Matrix::solve(Number& n)
{

	if (this->column_count() > this->row_count() + 1)
	{
		string message = "matrix is not solvable";
		
		cout << message << endl; throw string(message);
	}
	
	Matrix M (*this);
	
	
	//  Solve the matrix for X using Gaussian elimination
	
	//  First reduce the system to echelon form
	
	if (n.equals(0))                 M = M .to_echelon_form();
	else if (n.trim().length() >  1) M = M .to_echelon_form(n);
	else if (n.trim().length() == 1) M = M .to_echelon_form(n.int_value());
	
	
	if (!M .is_echelon_form())
	{
		string message = "system is not reducible to echelon form";
		
		for (int i = 0; i < M.column_count() -1; i++)
		
		    if (M.matrix[i][i]->equals(0))
		
			cout << "diagonal zero at row " << i;
		
		throw string(message);
	}
	
	
	//  Verify that the diagonal elements are non-zero
	
	Number zero(0); if (M .get_precision() != 0)
	
	    zero = zero .set_precision(M .get_precision());
	
	for (int i = 0; i < M .column_count() -1; i++)
	
	    if (M.get(i, i) .equals(zero))
	
		cout << "diagonal zero at row " << i;
	
	
	//  Reduce the system to row canonical form
	
	if (n.equals(0))  M = M .to_row_canonical_form();
	else              M = M .to_row_canonical_form(n);
	
	
	//  Verify that the system is solved
	//
	//  The matrix is solved if the left matrix of A|B is the unit matrix
	
	
	//  Copy the solution from the last column
	
	int length = M.column_count() -1;
	
	vector<Number*> column = M .get_column(length);
	
	vector<Number*> solution (length);
	
	for (int i = 0; i < length; i++)
	
	    solution[i] = column[i];
	
	if (!M .is_row_canonical_form())
	{
		//  if there are fewer equations than there are variables,
		//  or if some of the equations are not linearly independent
		//  then the system is unsolvable
		
		int digits = ((int) n .bit_count() + 3) % 4, radix = 10;
		
		cout << "System cannot be reduced to row canonical form";
		
		return solution;
	}
	
	//  Verify that A X == B
	
	cout << "Verifying A X == B";
	
	Matrix A = this->get(0, 0, length, length);
	
	//  Convert the solution to a column matrix
	
	Matrix X (solution);  X = X.transpose();
	
	//  Truncate the b column so that it equals the
	//  number of columns - 1 instead of the number of rows
	
	vector<Number*> b (length);
	
	for (int i = 0; i < length; i++)
	
	    b[i] = & this->get(i, this->column_count() -1);
	
	//  Convert the b vector to a column matrix
	
	Matrix B = Matrix(b) .transpose();
	
	Matrix A_X = A .multiply(X);
	
	if (!n.equals(0))
	{
		A_X = A_X .mod(n) .add(n) .mod(n);
		B   = B   .mod(n) .add(n) .mod(n);
	}
	
	if (!A_X .equals(B))
	{
		throw string("error");
	}
	
	return solution;
}


Matrix& Matrix::square()
{
	//  squares a matrix
	
	return this->multiply(*this);
}

Matrix& Matrix::subtract(int n)
{
	//  subtracts an int from a matrix
	
	Number number(n);
	
	return this->subtract(number);
}


Matrix& Matrix::subtract(Number& n)
{
	//  subtracts a number from a matrix
	
	Matrix matrix(*this);
	
	for (int i = 0; i < matrix.   row_count(); i++)
	for (int j = 0; j < matrix.column_count(); j++)
	
	    matrix.matrix[i][j] = &
	    matrix.matrix[i][j] ->subtract(n);
	
	return set_result(matrix);
}


Matrix& Matrix::subtract(Matrix& matrix)
{
	//  subtracts two matrices
	
	Matrix matrix1(*this);
	
	for (int i = 0; i < matrix1.row_count(); i++)
	for (int j = 0; j < matrix1.column_count(); j++)
	
	    matrix1.matrix[i][j] = &
	    matrix1.matrix[i][j]
	
		->subtract(*matrix.matrix[i][j]);
	
	return set_result(matrix1);
}


Matrix& Matrix::to_echelon_form()
{

	Matrix matrix(*this);
	
	int p = matrix.get_precision();
	
	//  Set a minimum precision to avoid dividing by zero
	
	if (p == 0)  p = 8;
	
	matrix = matrix .set_precision(p);
	
	int rows = matrix .row_count();
	
	cout << endl << matrix.to_matrix_string(16) << endl << endl;
	
	for (int r = 0; r < rows -1; r++)
	{
		//  Find the row i and column j of the first non-zero element
		
		bool bool1 = false;  int i = 0, j = 0;
		
		for (j = r; j < matrix.column_count(); j++)
		{
			for (i = r; i < matrix.row_count(); i++)
			
			    if (!matrix.matrix[i][j] ->equals(0.0))
			
				{ bool1 = true;  break; }
			
			if (bool1 == true)  break;
		}
		
		if (j >= matrix.column_count()) continue;
		
		//  Swap rows i and r if element [r][j] equals zero
		
		if (i != r)
		{
			vector<Number*> tempi = matrix.matrix[i];
			vector<Number*> tempr = matrix.matrix[r];
			
			matrix.matrix[i] = tempr;
			matrix.matrix[r] = tempi;
		}
		
		
		//  Put zeros below the pivot for each row i below r
		
		for (i = r + 1; i < rows; i++)
		{
			//  For each row i compute the multiplier
			//
			//  m = - [i][j] / [r][j]
			
			//  Possible divide by zero error if precision equals
			//  zero depending on the implementation of the inverter.
			//  This is avoided by setting a minimum precision.
			
			Number m = matrix.matrix[i][j] ->negate() .multiply(
			           matrix.matrix[r][j] ->inverse() );
			
			//  Multiply row r by m, and add to row i
			
			for (int k = j; k < matrix.column_count(); k++)
			
			    matrix.matrix[i][k] = & matrix.matrix[i][k]
			
				->add( matrix.matrix[r][k]->multiply(m) );
		}
	}
	
	matrix = matrix.set_precision(p);
	
	cout << endl << matrix.to_matrix_string(16) << endl << endl;
	
	return set_result(matrix);
}


Matrix& Matrix::to_echelon_form(int n)
{

	//  reduces a matrix to (upper) echelon form
	
	
	Matrix matrix(*this);
	
	matrix = matrix.mod(n);
	
	int rows = matrix.row_count();
	
	for (int r = 0; r < rows -1; r++)
	{
		//  Reduce the augmented matrix [ A | B ] to echelon form
		
		//  Find the row i and column j of the first non-zero element
		
		bool bool1 = false;  int i = 0, j = 0;
		
		for (j = r; j < matrix.column_count(); j++)
		{
			for (i = r; i < matrix.row_count(); i++)
			{
				if ((matrix.matrix[i][j]->length() == 1)
				 && (matrix.matrix[i][j]->to_int_vector()[0] == 0))
				
				    continue;
				
				if (!matrix.matrix[i][j]->equals(0))
				
				    { bool1 = true;  break; }
			}
			
			if (bool1 == true)  break;
		}
		
		if (j >= matrix.column_count()) continue;
		
		//  Swap rows i and r if element [r][j] equals zero
		
		if (i != r)
		{
			vector<Number*> tempi = matrix.matrix[i];
			vector<Number*> tempr = matrix.matrix[r];
			
			matrix.matrix[i] = tempr;
			matrix.matrix[r] = tempi;
		}
		
		//  Put zeros below the pivot for each row i below r
		
		for (i = r + 1; i < rows; i++)
		{
			//  For each row i compute the multipliers
			//
			//  m1 = - [i][j],  m2 = [r][j]
			//
			//  instead of  m = - [i][j] / [r][j]
			//
			//  to avoid the modular inversion.
			//
			//  [i][j] m2 + [r][j] m1 == [i][j][r][j] - [i][j][r][j] == 0
			//
			//  This method of elimination can only be used for modular systems
			//  because the coefficient size grows polynomially or quadratically
			
			Number m1 =  matrix.matrix[i][j] ->negate(n);
			Number m2 = *matrix.matrix[r][j];
			
			int m1int = m1.int_value(), m2int = m2.int_value();
			
			//  Multiply row r by m1, and add to row i multiplied by m2
			
			for (int k = j; k < matrix.column_count(); k++)
			
			    matrix.matrix[i][k] = new Number( ( 1L * matrix.matrix[i][k]->int_value() * (m2int)
			
				+ ( 1L * matrix.matrix[r][k]->int_value() * (m1int) ) ) % n );
		}
	}
	
	matrix = matrix .mod(n) .add(n) .mod(n);
	
	return set_result(matrix);
}


Matrix& Matrix::to_echelon_form(Number& n)
{

	//  reduces a matrix to echelon form
	
	//  The matrix can be square or rectangular (as in the determinant
	//  of a square matrix (n, n), the inverse of a rectangular matrix
	//  A X == I (n, 2 n), or the linear system A X == B  (n, n + 1)).
	//
	//  A matrix in echelon or triangular form is a matrix in which
	//  the first non-zero element of each row is to the right of the
	//  first non-zero element in the preceding row.
	
	
	//  Example of an (n, n + 1) matrix in echelon form
	//
	//  |  3  1  4  1  5  9  2  6  5  3 | 5 |
	//  |     8  9  7  9  3  2  3  8  4 | 6 |
	//  |        2  6  4  3  3  8  3  2 | 7 |
	//  |           9  5  0  2  8  8  4 | 1 |
	//  |              9  7  1  6  9  3 | 9 |
	//  |                 9  3  7  5  1 | 0 |
	//  |                    5  8  2  0 | 9 |
	//  |                       7  4  9 | 4 |
	//  |                          4  5 | 9 |
	//  |                             2 | 3 |
	
	
	
	//  Example of an (n, 2 n) matrix in echelon form
	//
	//  | 3  1  4  1  5  9 | 2  6  5  3  5  8 |
	//  |    9  7  9  3  2 | 3  8  4  6  2  6 |
	//  |       4  3  3  8 | 3  2  7  9  5  0 |
	//  |          2  8  8 | 4  1  9  7  1  6 |
	//  |             9  3 | 9  9  3  7  5  1 |
	//  |                5 | 8  2  0  9  7  4 |
	
	
	
	//  If the modulus is an int, use the int method
	
	if (n.length() == 1) return
	
	    to_echelon_form(n.int_value());
	
	
	
	//  Forward Elimination
	//
	//  Reduce the matrix to (upper) echelon form
	//
	//  (put zeros below each pivot)
	
	
	//  Pre-compute the inverse for fast modular reduction
	
	int digits = (int) n.bit_count() / 4;
	
	Number invn = n .set_precision(16 + digits) .inverse();
	
	
	Matrix matrix(*this);
	
	matrix = matrix .mod(n);
	
	int rows = matrix.row_count();
	
	for (int r = 0; r < rows - 1; r++)
	{
		if (matrix.row_count() > 256) cout << r << " ";
		
		//  Reduce the augmented matrix [ A | B ] to echelon form
		
		//  Find the row i and column j of the first non-zero element
		
		bool bool1 = false;  int i = 0, j = 0;
		
		for (j = r; j < matrix.column_count(); j++)
		{
			for (i = r; i < matrix.row_count(); i++)
			{
				if ((matrix.matrix[i][j]->length() == 1)
				 && (matrix.matrix[i][j]->to_int_vector()[0] == 0))
				
				    continue;
				
				if (!matrix.matrix[i][j]->equals(0))
				
				    { bool1 = true;  break; }
			}
			
			if (bool1 == true)  break;
		}
		
		if (j >= matrix.column_count()) continue;
		
		
		//  Swap rows i and r if element [r][j] equals zero
		
		if (i != r)
		{
			vector<Number*> tempi = matrix.matrix[i];
			vector<Number*> tempr = matrix.matrix[r];
			
			matrix.matrix[i] = tempr;
			matrix.matrix[r] = tempi;
		}
		
		
		//  Put zeros below the pivot for each row i below r
		
		for (i = r + 1; i < rows; i++)
		{
			//  For each row i compute the multipliers
			//
			//  m1 = - [i][j],  m2 = [r][j]
			//
			//  instead of  m = - [i][j] / [r][j]
			//
			//  to avoid the modular inversion.
			//
			//  [i][j] m2 + [r][j] m1 == [i][j][r][j] - [i][j][r][j] == 0
			//
			//  This method of elimination can only be used for modular systems
			//  because the coefficient size grows polynomially or quadratically
			
			Number m1 =  matrix.matrix[i][j] ->negate(n);
			Number m2 = *matrix.matrix[r][j];
			
			
			//  Multiply row r by m1, and add to row i multiplied by m2
			
			for (int k = j; k < matrix.column_count(); k++)
			
			    matrix.matrix[i][k] = & matrix.matrix[i][k]->multiply(m2)
			
				.add( matrix.matrix[r][k]->multiply(m1) ) .mod(n, invn);
		}
	}
	
	return set_result(matrix);
}



string& Matrix::to_integer_string(int digits, int radix)
{

	//  converts a matrix to an integer string by concatenating
	//  the elements. This method is used for cryptography to
	//  convert a public key (matrix) to a public key string.
	
	//  The digits variable is the minimum number of digits.
	//  The left side will be padded with zeros if necessary.
	
	//  This method is the inverse of the constructor Matrix(
	//  String str, int radix, int rows, int columns)
	
	
	Matrix matrix(*this);
	
	int rows    = matrix.row_count();
	int columns = matrix.column_count();
	
	string y[rows][columns];
	
	for (int i = 0; i < rows;    i++)
	for (int j = 0; j < columns; j++)
	{
		y[i][j] = matrix.matrix[i][j]
		
		    ->to_string(digits, radix);
		
		int pos = 0;
		
		while ((pos = y[i][j] .find_first_of(" ")) != string::npos)
		
		    y[i][j] = y[i][j] .erase(pos, 1);
		
		//  the toString method is guaranteed to prepend zeros
		//  if (radix == 16) unless the implementation changes
		
		//  Replace spaces with zeros
		
		if (digits == 0) continue;
		
		while (y[i][j].length() < digits)
		
		       y[i][j] = "0" + y[i][j];
	}
	
	string numberstr ("");
	
	for (int i = 0; i < rows; i++)
	for (int j = 0; j < columns; j++)
	
	    numberstr += y[i][j];
	
	return *new string(numberstr);
}


string& Matrix::to_matrix_string()
{
	return to_matrix_string(10);
}

string& Matrix::to_matrix_string(int radix)
{
	return to_matrix_string(0, radix);
}

string& Matrix::to_matrix_string(int digits, int radix)
{
	//  converts a matrix to a string that can
	//  be displayed as a rectangular matrix
	//
	//  | a11  a12  a13  a14  ... |
	//  | a21  a22  a23  a24  ... |
	//  | a31  a32  a33  a34  ... |
	//  | a41  a42  a43  a44  ... |
	//  | ...  ...  ...  ...  ... |
	//
	//  instead of the to_string() method which converts the matrix to
	//
	//  { { a11, a12, a13, ... }, { a21, a22, a23, ... }, ... }.
	
	
	int mindigits = 0, digitsize = 0;
	
	int rows = this->row_count();
	
	int columns = this->column_count();
	
	for (int i = 0; i < rows; i++)
	for (int j = 0; j < columns; j++)
	{
		if (this->matrix[i][j] ->equals(0))
		
		    { if (4 > mindigits) mindigits = 4;  continue; }
		
		string str = this->matrix[i][j] ->to_string(radix);
		
		int pos = 0;
		
		while ((pos = str.find_first_of(" ")) != string::npos)
		
		    str = str .erase(pos, 1);
		
		digitsize = str.length() + ((this->matrix[i][j]->signum() == -1) ? 1 : 0);
		
		if (digitsize > mindigits)  mindigits = digitsize;
	}
	
	string sb(""), sbdigit;
	
	for (int i = 0; i < rows; i++)
	{
		sb.append("| ");
		
		for (int j = 0; j < columns; j++)
		{
			string numberstr;
			
			//  if (this->matrix[i][j] != null)
			
			numberstr = this->matrix[i][j] ->to_string(digits, radix);
			
			int pos = 0;
			
			while ((pos = numberstr.find_first_of(" ")) != string::npos)
			
			    numberstr = numberstr .erase(pos, 1);
			
			sbdigit = string(numberstr);
			
			while (sbdigit.length() < mindigits)
			
			    sbdigit.insert(0, " ");
			
			numberstr = string(sbdigit);
			
			sbdigit = string("");
			
			sbdigit .append(numberstr);
			
			
			//  Append digit and spaces
			
			sb.append(sbdigit);
			
			int spaces = mindigits - (!numberstr.empty() ?
			
			    numberstr.length() : 4);
			
			if (j < columns -1) sb.append("  ");
			
			for (int k = 0; k < spaces; k++)  sb.append(" ");
		}
		
		sb .append(" |");
		
		if (i < rows -1)  sb.append("\n");
	}
	
	sb .append("\n");
	
	return *new string(sb);
}


Matrix& Matrix::to_row_canonical_form()
{
	Number zero(0);
	
	return this->to_row_canonical_form(zero);
}


Matrix& Matrix::to_row_canonical_form(Number& n)
{
	if (!is_echelon_form())
	{
		string message = "matrix is not in echelon form";
		
		cout << message << endl; throw string(message);
	}
	
	//  Determine which method to use
	
	//  If the matrix is augmented with a column, then use method one
	//  If the matrix is augmented with a matrix, then use method two
	
	int method = (column_count() <= row_count() + 1) ? 1 : 2;
	
	if   (method == 1) return to_row_canonical_form1(n);
	else               return to_row_canonical_form2(n);
}


Matrix& Matrix::to_row_canonical_form1(Number& n)
{

	//  reduces a matrix in echelon form to row
	//  canonical form using back substitution
	
	cout << "to_row_canonical_form1 method" << endl;
	
	Matrix matrix (*this);
	
	if (!n.equals(0)) matrix = matrix .mod(n);
	
	int rows    = column_count() -1;
	int columns = column_count();
	
	int p = matrix .get_precision();
	
	vector<Number*> x (rows),  B (rows);
	
	vector<Number*> temp = matrix.get_column(matrix.column_count() -1);
	
	for (int i = 0; i < rows; i++)  B[i] = temp[i];
	for (int i = 0; i < rows; i++)  x[i] = new Number(0);
	
	for (int i = rows - 1; i >= 0; i--)
	{
		Number sum (0);
		
		for (int j = i + 1; j < rows; j++)
		{
			sum = sum .add(matrix.matrix[i][j]->multiply(*x[j]));
			
			if (!n.equals(0))  sum = sum .mod(n);
		}
		
		Number B1 = B[i] ->subtract(sum);
		
		if (!n.equals(0)) B1 = B1 .mod(n) .add(n) .mod(n);
		
		try
		{	Number inverse;
			
			if (matrix.matrix[i][i]->equals(0))
			{
				string message = "to_row_canonical_formMethod1 error";
				
				cout << message << endl; throw string(message);
			}
			
			try
			{	if (!n.equals(0)) inverse = matrix.matrix[i][i] ->mod_inverse(n);
				 else             inverse = matrix.matrix[i][i] ->    inverse();
			}
			
			catch (...)
			{
				string message = "to_row_canonical_form1 inversion error";
				
				cout << message << endl; throw string(message);
			}
			
			if (!n.equals(0))  x[i] = & B1 .multiply(inverse) .mod(n);
			else               x[i] = & B1 .multiply(inverse);
			
			
			if (!n.equals(0))  x[i] = & x[i] ->mod(n) .add(n) .mod(n);
		}
		
		catch (...)
		{
			//  If an element is non-invertible modulo n, then the
			//  matrix is not reducible to row canonical form
			
			string message = "matrix is not reducible to canonical form modulo n";
			
			cout << message << endl; throw string(message);
		}
	}
	
	//  Restore the precision of x because the inverse method may increase the precision
	
	for (int i = 0; i < rows; i++) x[i] = & x[i]->set_precision(p);
	
	//  Return the augmented identity matrix
	
	return identity_matrix(rows) .augment(x);
}




Matrix& Matrix::to_row_canonical_form2(int n)
{

	//  reduces a matrix in echelon form to row canonical
	//  form using backward elimination, not back substitution
	
	Matrix matrix(*this);
	
	matrix = matrix .trim() .mod(n);
	
	if (!matrix.is_echelon_form())
	{
		string message = "matrix is not in echelon form";
		
		cout << message << endl; throw string(message);
	}
	
	if (matrix.column_count() != matrix.row_count() + 1)
	{
		string message = "columns != rows + 1";
		
		cout << message << endl; throw string (message); // wrong method
	}
	
	Number zero (0);
	
	Matrix temp = matrix .get(0, 0,
	
	     matrix.row_count(), matrix.column_count());
	
	
	//  Divide each row so that the fist nonzero element equals one
	
	for (int r = matrix.row_count() -1, j = 0; r >= 0; r--)
	{
		//  Find the column of the first non-zero element
		
		for (j = 0; j < matrix.column_count() -1; j++)
		
		    if (!matrix.matrix[r][j]->equals(zero))  break;
		
		if (j >= matrix.column_count() -1)  continue;
		
		if (matrix.matrix[r][j]->mod(n) .equals(zero)) continue;
		
		
		//  Compute the multiplier m = 1 / pivot
		
		Number inv;
		
		try { inv = matrix.matrix[r][j]->mod_inverse(n); }
		
		catch (...)
		{
			string message = "to_row_canonical_form2 inversion error";
			
			cout << message << endl; throw string(message);
		}
		
		
		//  Multiply row r by the multiplier so the pivot equals one
		
		for (int k = 0; k < matrix.column_count(); k++)
		{
			matrix.matrix[r][k] = & matrix.matrix[r][k]->multiply(inv);
			
			matrix.matrix[r][k] = & matrix.matrix[r][k]->mod(n);
		}
		
		
		//  Put zeros above the pivot for each row above r
		
		for (int i = r - 1; i >= 0; i--)
		{
			//  For each row i compute the multipliers
			//
			//  m1 = - [i][j],   m2 = [r][j]
			//
			//  instead of  m = -[i][j] / [r][j]
			//
			//  to avoid the modular inversion.
			//
			//  [i][j] m2 + [r][j] m1 == [i][j][r][j] - [i][j][r][j] == 0
			
			
			int m1 = matrix.matrix[i][j]->int_value();
			int m2 = matrix.matrix[r][j]->int_value();
			
			m1 = n - m1;
			
			
			//  Multiply row r by m1, and add to row i multiplied by m2
			
			for (int k = j; k < matrix.column_count(); k++)
			{
				matrix.matrix[i][k] = new Number(
				
				     1L * matrix.matrix[i][k]->int_value() * m2
				 + ( 1L * matrix.matrix[r][k]->int_value() * m1 ) );
				
				matrix.matrix[i][k] = & matrix.matrix[i][k]->mod(n);
			}
		}
	}
	
	matrix = matrix .mod(n) .add(n) .mod(n);
	
	return set_result(matrix);
}



Matrix& Matrix::to_row_canonical_form2(Number& n)
{

	//  reduces a matrix in echelon form to row canonical
	//  form using backward elimination, not back substitution
	
	if (!is_echelon_form())
	{
		string message = "matrix is not in echelon form";
		
		cout << message << endl; throw string(message);
	}
	
	Matrix matrix(*this);
	
	if (!n.equals(0))  matrix = matrix.mod(n);
	
	else // if (!n.equals(0))
	{
		int precision = matrix.get_precision();
		
		if (precision == 0) precision = 8;
		
		matrix = matrix.set_precision(precision);
	}
	
	Number zero(0);
	
	int precision = matrix.get_precision();
	
	if (precision > 0)
	
	    zero = zero .set_precision(precision);
	
	
	//  Divide each row so that the fist non-zero element equals one
	
	for (int r = matrix.row_count() -1, j = 0; r >= 0; r--)
	{
		//  Find the column of the first non-zero element
		
		for (j = 0; j < matrix.column_count() -1; j++)
		
		    if (!matrix.matrix[r][j]->equals(zero))  break;
		
		if (j >= matrix.column_count() -1)  continue;
		
		
		if ((!n.equals(0)) && matrix.matrix[r][j]->mod(n) .equals(zero)) continue;
		if (( n.equals(0)) && matrix.matrix[r][j]        ->equals(zero)) continue;
		
		
		//  Compute the multiplier m = 1 / pivot
		
		Number inv;
		
		try
		{	if (!n.equals(0)) inv = matrix.matrix[r][j]->mod_inverse(n);
			 if ( n.equals(0)) inv = matrix.matrix[r][j]    ->inverse();
		}
		
		catch (...)
		{
			string message = "to_row_canonical_form2 inversion error";
			
			cout << message << endl; throw string(message);
		}
		
		
		//  Multiply row r by the multiplier so the pivot equals one
		
		for (int k = 0; k < matrix.column_count(); k++)
		{
			matrix.matrix[r][k] = & matrix.matrix[r][k]->multiply(inv);
			
			if (!n.equals(0)) matrix.matrix[r][k] = & matrix.matrix[r][k]->mod(n);
		}
		
		
		//  Put zeros above the pivot for each row above r
		
		for (int i = r-1; i >= 0; i--)
		{
			//  For each row i compute the multipliers
			//
			//  m1 = - [i][j],   m2 = [r][j]
			//
			//  instead of  m = -[i][j] / [r][j]
			//
			//  to avoid the modular inversion.
			//
			//  [i][j] m2 + [r][j] m1 == [i][j][r][j] - [i][j][r][j] == 0
			
			
			Number m1, m2;
			
			m1 = *matrix.matrix[i][j];
			m2 = *matrix.matrix[r][j];
			
			if (!n.equals(0)) m1 = m1 .negate(n);
			if ( n.equals(0)) m1 = m1 .negate();
			
			
			//  Multiply row r by m1, and add to row i multiplied by m2
			
			for (int k = j; k < matrix.column_count(); k++)
			{
				matrix.matrix[i][k] = & matrix.matrix[i][k]->multiply(m2)
				                  .add( matrix.matrix[r][k]->multiply(m1) );
				
				if (!n.equals(0)) matrix.matrix[i][k] = & matrix.matrix[i][k]->mod(n);
			}
		}
	}
	
	if (!n.equals(0)) matrix = matrix
	
	    .mod(n) .add(n) .mod(n);
	
	return set_result(matrix);
}


string& Matrix::to_string()
{
	return to_string(10);
}

string& Matrix::to_string(int radix)
{
	return to_string(0, radix);
}

string& Matrix::to_string(int digits, int radix)
{
	//  converts the matrix to string
	
	string str("");  str += "{";
	
	int rows = this->row_count();
	
	int columns = this->column_count();
	
	for (int i = 0; i < rows; i++)
	{
		str += " { ";
		
		for (int j = 0; j < columns; j++)
		{
			str += this->matrix[i][j] ->to_string(digits, radix);
			
			if (j < columns -1)  str += ", ";
		}
		
		str += " }";
		
		if (i < rows -1) str += ",";
	}
	
	str += " }";
	
	return *new string(str);
}



Number& Matrix::trace()
{
	//  returns the trace of a matrix
	//
	//  The trace is the sum of the diagonal components
	
	if (row_count() != column_count())
	{
		string message = "non-square matrix";
		
		cout << message << endl; throw string(message);
	}
	
	Number trace(0);
	
	for (int i = 0; i < row_count(); i++)
	
	    trace = trace .add(*this->matrix[i][i]);
	
	return *new Number(trace);
}


Matrix& Matrix::transpose()
{
	//  swaps the rows and columns of a matrix
	
	int rows    = this->row_count();
	int columns = this->column_count();
	
	Matrix matrix (columns, rows);
	
	for (int i = 0; i < rows;    i++)
	for (int j = 0; j < columns; j++)
	
	    matrix.matrix[j][i] = this->matrix[i][j];
	
	return set_result(matrix);
}


Matrix& Matrix::trim()
{
	//  removes rows of zeros from the bottom
	
	Matrix matrix(*this);
	
	for (int i = matrix.row_count() -1; i >= 0; i--)
	{
		vector<Number*> row = matrix.get_row(i);
		
		bool zero = true;
		
		for (int j = 0; j < column_count(); j++)
		
		    if (!row[j]->equals(0)) { zero = false;  break; }
		
		if (zero) matrix = matrix .delete_row(i);
		
		else break;
	}
	
	return set_result(matrix);
}













//  Copy control for the Number class
//
//  Define the copy constructor and assignment operator =


//  Note that the compiler has the option to rewrite an assignment as a
//  constructor call which uses the copy constructor instead of the
//  assignment operator. The gnu compiler calls the assignment operator
//  method if a variable is being reassigned such as a = a .function(),
//  otherwise it calls the copy constructor. For example, Number b = a;
//  is rewritten as Number b (a); which means that the copy constructor
//  is called instead of the assignment method.



//  Copy constructor

Number::Number(Number& original)
{
	//  cout << "number copy constructor" << endl;
	
	this->vec_int   = original.vec_int;
	this->intpoint  = original.intpoint;
	this->precision = original.precision;
	this->sign      = original.sign;
}



//  Assignment operator

Number& Number::operator = (Number& original)
{
	//  cout << "number assignment operator = " << endl;
	
	this->vec_int   = original.vec_int;
	this->intpoint  = original.intpoint;
	this->precision = original.precision;
	this->sign      = original.sign;
	
	return *this;
}












//  Number destructor
//
//
//  The Number destructor is called recursively for equations such as
//
//  A = A .function() .function() .function() .function() ...
//
//  because it creates pointers to pointers to pointers ... to Number,
//  or a linked list of pointers to Number.
//
//  As each function is evaluated, it returns a new anonymous object
//  which is used to evaluate the next function, which returns another
//  new anonymous object which is then used to evaluate the next function,
//  ... until the last new object is assigned to the variable.
//
//  Note that the Number class is not thread safe. If a Number is used
//  in the capture clause of a lambda function, the function has to make
//  a copy of the original Number using a constructor such as Number p1 (p)
//  and use the copy to prevent multiple threads from deleting the same
//  pointers which will cause a double free or memory corruption error.



Number::~Number()
{
	//  cout << "number destructor" << endl;
	
	//  Delete the pointers to pointers to Number
	//     and then delete the pointers to Number
	
	for (int i = 0; i < this->result.size(); i++)
	{
		if (this->result[i] != nullptr)
		{
			for (int j = 0; j < this->result[i]->result.size(); j++)
			{
				if (this->result[i]->result[j] != nullptr)
				{
					//  cout << "deleting pointer to pointer";
					
					delete this->result[i]->result[j];
			        	       this->result[i]->result[j] = nullptr;
				}
			}
			
			delete this->result[i];
			       this->result[i] = nullptr;
		}
	}
}




//  Number constructors


Number::Number()
{

	vector<int> vec_int { 0 };
	
	this->vec_int   = vec_int;
	this->intpoint  = 0;
	this->precision = 0;
	this->sign      = '+';
}


Number::Number(int val)
{
	char sign;
	
	if (val < 0)
	{
		//  two's complement
		
		val = ~val + 1;
		
		sign = '-';
	}
	
	else { sign = '+'; }
	
	vector<int> vec_int { val };
	
	this->vec_int   = vec_int;
	this->intpoint  = 0;
	this->precision = 0;
	this->sign      = sign;
}


Number::Number(long val)
{
	char sign;
	
	if (val < 0)
	{
		//  two's complement
		
		val = ~val + 1;
		
		sign = '-';
	}
	
	else { sign = '+'; }
	
	
	vector<int> vec_int(2);
	
	vec_int[0] = val >> 32;
	vec_int[1] = val >>  0;
	
	this->vec_int   = vec_int;
	this->intpoint  = 0;
	this->precision = 0;
	this->sign      = sign;
}


Number::Number(double d)
{
	//  Convert the double value to a Number
	
	Number *number = double_to_number(d);
	
	//  Make a copy of the number int vector
	
	int size = number->vec_int.size();
	
	vector<int> vec_int (size);
	
	for (int i = 0; i < size; i++)
	
	    vec_int[i] = number->vec_int[i];
	
	this->vec_int   = vec_int;
	this->intpoint  = number->intpoint;
	this->precision = number->precision;
	this->sign      = number->sign;
	
	delete number;
}


Number::Number(string digits)
{

	//  radix == 10
	
	//  Convert the string of digits to a Number
	
	Number *number = string_to_number(digits, 10);
	
	int size = number->vec_int.size();
	
	vector<int> vec_int (size);
	
	for (int i = 0; i < size; i++)
	
	    vec_int[i] = number->vec_int[i];
	
	this->vec_int   = vec_int;
	this->intpoint  = number->intpoint;
	this->precision = number->precision;
	this->sign      = number->sign;
	
	delete number;
}


Number::Number(string digits, int radix)
{
	Number *number = string_to_number(digits, radix);
	
	int size = number->vec_int.size();
	
	vector<int> vec_int (size);
	
	for (int i = 0; i < size; i++)
	
	    vec_int[i] = number->vec_int[i];
	
	this->vec_int   = vec_int;
	this->intpoint  = number->intpoint;
	this->precision = number->precision;
	this->sign      = number->sign;
	
	delete number;
}


Number::Number(vector<int> vec_int)
{
	int size = vec_int.size();
	
	vector<int> vec_int1 (size);
	
	for (int i = 0; i < size; i++)
	
	    vec_int1[i] = vec_int[i];
	
	this->vec_int   = vec_int1;
	this->intpoint  = 0;
	this->precision = 0;
	this->sign      = '+';
}


Number::Number(vector<int> vec_int, bool signed1)
{
	int size = vec_int.size();
	
	vector<int> vec_int1 (size);
	
	char sign;
	
	if (signed1 == false)
	{
		for (int i = 0; i < size; i++)
		
		    vec_int1[i] = vec_int[i];
		
		sign = '+';
	}
	
	else // if (signed1) // two's complement
	{
		//  Complement the bits and add 1
		
		if ((size != 0) && (vec_int[0] < 0))
		{
			vec_int1 = Math::twos_complement(vec_int);
			
			sign = '-';
		}
		
		else
		{	for (int i = 0; i < size; i++)
			
			    vec_int1[i] = vec_int[i];
			
			sign = '+';
		}
	}
	
	this->vec_int   = vec_int1;
	this->intpoint  = 0;
	this->precision = 0;
	this->sign      = sign;
}



Number::Number(Number& real, Number& imag)
{
	//  This constructor creates a complex number
	
	if (real.is_complex() || imag.is_complex())
	{
		string message = "real and imag components are imaginary";
		
		cout << message << endl; throw string(message);
	}
	
	//  Real component
	
	int size = real.vec_int.size();
	
	vector<int> vec_int (size);
	
	for (int i = 0; i < size; i++)
	
	    vec_int[i] = real.vec_int[i];
	
	this->vec_int   = vec_int;
	this->intpoint  = real.intpoint;
	this->precision = real.precision;
	this->sign      = real.sign;
	
	
	//  Imag component
	
	int size1 = imag.vec_int.size();
	
	vector<int> vec_int1 (size);
	
	for (int i = 0; i < size1; i++)
	
	    vec_int1[i] = imag.vec_int[i];
	
	this->vec_int1    = vec_int1;
	this->intpoint1  = imag.intpoint;
	this->precision1 = imag.precision;
	this->sign1      = imag.sign;
}




Number* Number::double_to_number(double d)
{
	string str = std::to_string(d);
	
	int index = -1;
	
	for (int i = 0; i < str.length(); i++)
	
	    if (isalpha(str[i])) // e index
	
		{ index = i; break; }
	
	int exp = 0;  // e + xx
	
	string substr = str.substr(index + 1);
	
	if (index != -1) exp = parse_int(substr);
	
	string significand = ((index != -1) ?
	
	    str.substr(0, index) : str);
	
	if (significand.find(".") != string::npos)
	
	    while (significand.length() < 15)
	
		significand += string("00");
	
	//  Set the number equal to the significand
	
	Number *number = string_to_number(significand, 10);
	
	//  Set the integer fraction point
	
	string base_10 ("10");
	
	int precision = number->precision;
	
	
	Number* temp = string_to_number(base_10, 10);
	
	temp = & temp->set_precision(precision) .pow((int) exp);
	
	number = & number->multiply(*temp);
	
	delete temp;
	
	return number;
}





//  Static methods used by the Number constructors


Number* Number::string_to_number(string& digits, int radix)
{
	char sign = '+';
	
	string str = digits;
	
	for (auto & c : str)
	
	    c = tolower(c);
	
	int pos = 0; // pos of the integer point
	int exp = 0; // exp of the significand
	
	//  Trim the string
	
	while ((pos = str.find(" ")) != string::npos)
	
	    str = str .erase(pos, 1);
	
	if (str.empty())
	{
		string message = "empty constructor string";
		
		cout << message << endl; throw string(message);
	}
	
	if (str[0] == '-')
	{
		str = str.substr(1, str.length() - 1);
		
		sign = '-';
	}
	
	else if (str[0] == '+')
	{
		str = str.substr(1, str.length() - 1);
		
		sign = '+';
	}
	
	else { sign = '+'; }
	
	
	if (!is_number_string(str, radix))
	{
		string message = string("number constructor string ") + str
		
		 + string(" is not a number string in radix ")
		 
		     + Number(radix).to_string();
		
		cout << message << endl; throw string(message);
	}
	
	//  If the string ends with a period, remove the period;
	//  else if the string contains a period, set the exponent
	//  to the string length - 1 - the position of the period.
	
	if (str.find_last_of(".") == str.length() -1)
	
	    str = str .substr(0, str.length() -1);
	
	else if (str.find(".") != string::npos)
	{
		int pos = str.find(".");
		
		exp = str.length() -1 - pos;
		
		str = str.erase(pos, 1);
	}
	
	int precision = exp;
	
	Number *number = new Number(0);
	
	if (radix == 16)
	{
		//  Convert from string to 4-bit hex chars (log 16)
		//  Convert from 4-bit chars (log 16) to 8-bit chars (log 256)
		//  Convert from 8-bit chars to 32-bit (ints)
		
		while ((str.length() % 8) != 0) str = "0" + str;
		
		int str_length = str.length();
		
		char *char_16_array = Convert::string_to_char_array(str);
		
		int char_16_length = str_length;
		
		char *char_256_array = Convert::char_16_array_to_char_256_array(
		
		    char_16_array, char_16_length);
		
		delete[] char_16_array;
		
		int char_256_length = char_16_length / 2;
		
		int *int_array = Convert::char_array_to_int_array(
		
		    char_256_array, char_256_length);
		
		delete[] char_256_array;
		
		int int_array_length = char_256_length / 4;
		
		//  delete and reassign
		
		number->vec_int = vector<int> (int_array_length);
		
		for (int i = 0; i < int_array_length; i++)
		
		    number->vec_int[i] = int_array[i];
		
		delete[] int_array;
		
		number->precision = precision;
		
		//  Set the integer fraction point
		
		if (exp > 0) //  divide by 16 ^ exp (w/o using a divider)
		{
			int div = exp / 8, mod = exp % 8;
			
			number->intpoint += div;
			
			if (mod != 0)
			{
				number = & number->shift_left(32, 32);
				
				number->intpoint += 1;
				
				number = & number->shift_right(4 * mod);
			}
		}
	}
	
	else // if (radix != 16)
	{
		//  Compute the integer value
		
		delete number;
		
		number = string_to_number1(str, radix);
		
		//  Compute the integer fraction point
		
		if (exp > 0)
		{
			number->precision = precision;
			
			number = & number ->divide(
			
			    Number(radix).pow((int) exp));
		}
	}
	
	number->sign = sign;
	
	return number;
}



Number* Number::string_to_number1(string &digits, int radix)
{

	//  This method bifurcates an n-digit integer string into two d and
	//  (n - d) length strings where d = (approx) n / 2, computes the
	//  left and right numbers separately, and then computes the number
	//  n = (left x radix ^ d) + right. This reduces the running time
	//  from quadratic to linear log or O(n^2) to O(n log(n)).
	
	
	int maxsize = 256;
	
	int pos = 0;
	
	string str = digits;
	
	while ((pos = str.find(" ")) != string::npos)
	
	    str = str .erase(pos, 1);
	
	
	//  Test before recursion
	
	if (str.length() < maxsize)
	{
		//  Convert a small string to a number
		
		Number number(0);
		
		for (int i = 0; i < str.length(); i++)
		{
			char c = str[str.length() - 1 - i];
			
			int digit = 0;
			
			if ((c >= '0') && (c <= '9'))
			
			     digit = c - '0';
			
			else if ((c >= 'a') && (c <= 'z'))
			
			    digit = c + 10 - 'a';
			
			number = number .add( Number(radix)
			
			    .pow((int) i) .multiply(digit) );
		}
		
		return new Number (number);
	}
	
	//  else if (str.length() >= maxsize)
	
	//  Convert a large string to a number
	
	int d = str.length() / 2;
	
	string lstring = str.substr(0, str.length() - d);
	
	string rstring = str.substr(str.length() - d,
	
	    str.length() - (str.length() - d));
	
	
	//  Start recursion
	
	return & string_to_number1(lstring, radix)
	
	    ->multiply( Number(radix) .pow((int) d) )
	
		.add( *string_to_number1(rstring, radix) );
}





Number& Number::set_result(Number& number)
{
	Number *result = new Number(number);
	
	this->result.push_back(result);
	
	return *result;
}




//  Number methods



Number& Number::abs()
{
	//  returns the absolute value
	
	if (this->is_complex())
	
	    throw string("IllegalArgumentException");
	
	Number number (*this);
	
	if (number.sign == '-')
	
	    number.sign = '+';
	
	return set_result(number);
}


Number& Number::add(int addend)
{
	Number number(addend);
	
	return this->add(number);
}

Number& Number::add(long addend)
{
	Number number(addend);
	
	return this->add(number);
}

Number& Number::add(double addend)
{
	Number number(addend);
	
	return this->add(number);
}

Number& Number::add(Number& addend)
{
	//  adds two numbers
	
	Number a (*this);
	Number b (addend);
	
	int precision = (a.precision > b.precision) ?
	                 a.precision : b.precision;
	
	a = a.set_precision(precision);
	b = b.set_precision(precision);
	
	if (!this->is_integer() || !addend.is_integer())
	{
		//  Align the integer / fraction points
		
		const int r = a.intpoint - b.intpoint;
		
		if      (r > 0) { b = b.shift_left( 32*r,  32*r);  b.intpoint += r; }
		else if (r < 0) { a = a.shift_left(-32*r, -32*r);  a.intpoint -= r; }
	}
	
	//  Expand the vectors by one int for the sign
	
	a = a .shift_left(32, 0);
	b = b .shift_left(32, 0);
	
	//  Make copies of the vectors
	
	vector<int> a_signed = a.vec_int;
	vector<int> b_signed = b.vec_int;
	
	Number a1 (a_signed, true); a1.intpoint = intpoint; a1.precision = precision;
	Number b1 (b_signed, true); b1.intpoint = intpoint; b1.precision = precision;
	
	
	//  Convert the unsigned int vectors to signed vectors if they are not already signed
	
	if (a.sign == '-')  a_signed = Math::twos_complement(a_signed);
	if (b.sign == '-')  b_signed = Math::twos_complement(b_signed);
	
	
	//  Add the two signed vectors
	
	vector<int> c_signed = Math::add(a_signed, b_signed);
	
	
	//  If the zeroth element < 0, then the vector is signed
	
	Number c (c_signed, true);
	
	c.intpoint = a.intpoint;
	c.precision = precision;
	
	c = c.trim();
	
	
	//  Return the result
	
	if (!this->is_complex() && !addend.is_complex())
	
	    return set_result(c);
	
	
	//  If the number is not real then use complex addition
	
	Number real1 = this->to_real();
	Number imag1 = this->to_imag();
	
	Number real2 = addend.to_real();
	Number imag2 = addend.to_imag();
	
	Number complex_sum
	
	    (real1.add(real2),
	     imag1.add(imag2));
	
	return set_result(complex_sum);
}


Number& Number::add_bit(long bit)
{
	//  adds a bit to a number
	
	//  The number is changed by this method
	//
	//  and by the set_bit and clear_bit methods
	
	int this_length = this->vec_int.size();
	
	if ((bit / 32 + 1) > this_length)
	{
		//  Expand the vector
		
		int vector_length = (int) (bit) / 32 + 1;
		
		vector<int> vec_int (vector_length);
		
		for (int i = 0; i < this_length; i++)
		
		    vec_int[vector_length - this_length] = this->vec_int[i];
		
		this->vec_int = vec_int;
	}
	
	Math::add_bit(this->vec_int, bit);
	
	return  *this;
}


Number& Number::arc_cos()
{
	//  computes the inverse cosine
	
	int p = this->precision;
	
	Number arc_cos = Number::pi(p) .divide(2)
	
	    .subtract(this->arc_sin());
	
	return set_result(arc_cos);
}


Number& Number::arc_cosh()
{
	//  computes the inverse hyperbolic cosine
	
	//  ...
}

Number& Number::arc_sin()
{

	//  computes the inverse sine
	
	//  ...
}

Number& Number::arc_sinh()
{
	//  computes the inverse hyperbolic sine
	
	//  ...
}

Number& Number::arc_tan()
{
	//  computes the inverse tangent
	
	//  ...
}


Number& Number::arc_tanh()
{
	//  computes the inverse hyperbolic tangent
	
	//  ...
}

Number& Number::and_(Number& val)
{
	//  ands two numbers
	
	vector<int> vec_int = Math::and_(this->vec_int, val.vec_int);
	
	int vec_size = this->vec_int.size() < val.vec_int.size() ?
	                   this->vec_int.size() : val.vec_int.size();
	
	Number number(vec_int, vec_size);
	
	number.intpoint  = this->intpoint;
	number.precision = this->precision;
	number.sign      = this->sign;
	
	return set_result(number);
}

long Number::bit_count()
{
	//  returns the highest set bit + 1
	
	return Math::bit_count(this->vec_int);
}

Number& Number::clear_bit(long bit)
{
	//  The number is changed by the clear_bit method
	//  and by the set_bit and add_bit methods
	
	Math::clear_bit(this->vec_int, bit);
	
	return *this;
}

int Number::compare(double val)
{
	Number number(val);
	
	return this->compare(number);
}

int Number::compare(int val)
{
	Number number(val);
	
	return this->compare(number);
}

int Number::compare(long val)
{
	Number number(val);
	
	return this->compare(number);
}

int Number::compare(Number& val)
{
	//  compares two numbers and returns 1, 0, or -1
	
	Number a = this->trim();
	
	Number b = val.trim();
	
	//  If signs are opposite, return 1 or -1
	
	if      ((a.sign != '-') && (b.sign == '-'))  return  1;
	else if ((a.sign == '-') && (b.sign != '-'))  return -1;
	
	
	//  If signs are both negative, swap the two numbers
	
	if ((a.sign == '-') && (b.sign != '-')) 
	{
		Number temp = Number(a) .lvalue();
		
		a = b.abs();
		
		b = temp.abs();
	}
	
	
	//  This method is not finished
	
	//    ...
	
	//    ...
	
	
	return 0;
}


Number& Number::cos()
{
	Number cos = cos_sin()[0];
	
	return set_result(cos);
}

Number *Number::cos_sin()
{
	//  returns a cosine and sine pair
	
	//  The series for cosine and sine are
	//
	//                n  2 n + 0
	//  cos(x) == (-1)  x       / (2n + 0)!
	//
	//                n  2 n + 1
	//  sin(x) == (-1)  x       / (2n + 1)!
	
	
	const int bits = 16;
	
	int p = (this->precision > 8) ?
	         this->precision + bits : 8 + bits;
	
	Number x = Number(*this) .set_precision(p);
	
	
	//  Divide x by 2^k
	
	int k = Math::bit_count(x.int_value()) + bits;
	
	x = x .divide( Number(2).pow(k)
	
	    .set_precision(x.precision + k) );
	
	Number cosx(0);
	Number sinx(0);
	
	Number inv_divisor = Number(1).set_precision(p);
	
	for (int n = 0; n < this->precision; n++)
	{
		//  cosx = cosx.add( x.pow( 2*n + 0 )
		//
		//      .divide( factorial( 2*n + 0 ) )
		//
		//  	  .multiply(((n % 2) == 0) ? 1 : -1) );
		
		//  (2n + 0)! == 0!, 2!, 4!, 6!, ... == 1, 2, 24, 720, ...
		
		if (n > 0)
		{
			inv_divisor = inv_divisor .divide(2*n - 0);
			inv_divisor = inv_divisor .divide(2*n - 1);
		}
		
		cosx = cosx .add(x.pow(2*n + 0)
		
		    .multiply(inv_divisor) .multiply((n % 2) == 0 ? 1 : -1));
	}
	
	inv_divisor = Number(1).set_precision(p);
	
	for (int n = 0; n < this->precision; n++)
	{
		//  sinx = sinx.add( x.pow( 2*n + 1 )
		//
		//      .divide( factorial( 2*n + 1 ) )
		//
		//        .multiply(((n % 2) == 0) ? 1 : -1) );
		
		//  (2n + 1)! == 1!, 3!, 5!, 7!, ... == 1, 6, 120, 5040, ...
		
		if (n > 0)
		{
			inv_divisor = inv_divisor .divide(2*n + 1);
			inv_divisor = inv_divisor .divide(2*n + 0);
		}
		
		sinx = sinx .add(x.pow(2*n + 1)
		
		    .multiply(inv_divisor) .multiply((n % 2) == 0 ? 1 : -1));
	}
	
	//  Multiply x by 2^k using the double-angle formulas
	//
	//  cos(2a) == cos^2(a) - sin^2(a)
	//  sin(2a) == 2 sin(a) cos(a)
	
	for (int i = 0; i < k; i++)
	{
		Number cosx1 = cosx.square() .subtract(sinx.square());
		
		sinx = sinx .multiply(cosx) .multiply(2);
		cosx = cosx1;
	}
	
	cosx = cosx .set_precision(this->precision);
	
	Number *cos_sin = new Number[2];
	
	cos_sin[0] = cosx;  cos_sin[1] = sinx;
	
	return cos_sin;
}


Number& Number::cosh()
{

	//  The hyperbolic cosine function is
	//  
	//                       u      -u
	//  cosh(u) =  1 / 2  ( e   +  e  )
	//
	//  The hyperbolic cosine is related to the cosine by
	//
	//  cos(u) == cosh(i u)  and  cosh(u) == cos(i u)
	
	
	this->precision = (this->precision > 0) ? this->precision : 16;
	
	Number cosh = e_(this->precision) .pow(*this) .subtract(
	
	              e_(this->precision) .pow(this->negate())) .divide(2);
	
	return set_result(cosh);
}

Number& Number::cube()
{
	//  cubes a number
	
	Number n (*this);
	
	Number cube = n .multiply(n) .multiply(n);
	
	return set_result(cube);
}

Number& Number::divide(int divisor)
{
	Number number(divisor);
	
	return this->divide(number);
}

Number& Number::divide(long divisor)
{
	Number number(divisor);
	
	return this->divide(number);
}

Number& Number::divide(double divisor)
{
	Number number(divisor);
	
	return this->divide(number);
}

Number& Number::divide(Number& divisor)
{
	//  divides one number by another number
	
	//  If the numbers are both integers, then the quotient will be an integer.
	//  If the dividend or the divisor is a real number, then the quotient will
	//  be a real number.
	
	//  Division can be done by inversion and multiplication using Newton's
	//  iteration for inversion (short division), or by partial quotients
	//  using a quadratic divider (long division).
	//
	//  The divide method determines which divider to call depending on the
	//  sizes and relative sizes of the dividend and divisor.
	//
	//  For computing the value of pi >= 1024 digits, the inverter is several
	//  times faster than the quadratic divider, but for computing public keys,
	//  the quadratic divider is a few times faster than the inverter.
	
	
	if (divisor.equals(0))
	{
		if (this->equals(0)) // 0 / 0 == 1
		
		    return set_result(Number(1).lvalue());
		
		string message = "divide by zero";
		
		cout << message << endl; throw string(message);
	}
	
	if (this->equals(0) || divisor.equals(1))
	
		return set_result(*this);
	
	
	if (divisor.is_complex())
	{
	
		//  Compute the inverse of a complex number a + i b
		//
		//	 1     (a - i b)       a - i b
		//  ==  _______  _________  ==  _________
		//
		//      a + i b  (a - i b)      a^2 + b^2
		
		
		//  denominator = a^2 + b^2 == real^2 + imag^2
		
		if (divisor.precision == 0) divisor.precision = 8;
		
		Number a = divisor.to_real();
		Number b = divisor.to_imag();
		
		Number real = a .divide( a.square() .add(b.square()) );
		Number imag = b .divide( a.square() .add(b.square()) );
		
		Number inverse(real, imag.negate());
		
		//  this / complex divisor == this * inv(complex divisor)
		
		Number quotient = this->multiply(inverse);
		
		
		//  Verify the complex inverse
		
		if (!quotient .multiply(divisor) .equals(*this))
		{
			string message = "complex division error";
			
			cout << message << endl; throw string(message);
		}
		
		quotient = quotient .trim();
		
		return set_result(quotient);
	}
	
	
	if (this->is_complex())
	{
		Number real = this->to_real().divide(divisor);
		Number imag = this->to_imag().divide(divisor);
		
		return set_result(Number(real, imag).lvalue());
	}
	
	
	Number a =    this->abs() .trim();
	Number b = divisor .abs() .trim();
	
	
	char sign = (this->sign != divisor.sign) ? '-' : '+';
	
	int precision = (a.precision >= b.precision) ?
			 a.precision  : b.precision;
	
	
	//  Integer division
	
	
	if ((a.precision == 0) && (b.precision == 0))
	{
	
		Number q;
		
		int a_bits = (int) (a.bit_count());
		int b_bits = (int) (b.bit_count());
		
		int method = -1;
		
		
		//  Three special cases for division
		//  
		//  size a <= size b,  size a >> size b,  size a ~ size b
		
		
		if (a_bits <= b_bits)
		{
			q = Number(0).lvalue();
			
			if (a.sign == '-') q = q .negate();
			
			method = 0;
		}
		
		else if (a_bits - b_bits < 4)
		{
			//  Subtract the divisor from the dividend
			//    or add the divisor to the dividend
			
			//  This code makes the gcd method a few times faster
			//  because the dividend and divisor are approx the same size
			
			q = Number(0).lvalue();
			
			Number temp = a .abs();
			
			while (!temp .less_than(b))
			{
				temp = temp .subtract(b);
				
				q = q .add(1);
			}
			
			if (a.sign == '-') q = q .negate();
			
			method = 1;
		}
		
		//	else if ((a.length > 1) && (b.length == 1))
		//	{
		//		//  this method is not finished
		//		
		//		q = a .divide_by_int(b.int_value()) .to_integer();
		//		
		//		method = 2;
		//	}
		
		else if ((a_bits < 512) && (b_bits < 512))
		{
			//  quadratic division (long division)
			
			q = a.quad_divide(b);
			
			method = 3;
		}
		
		else
		{	//  Newton's Iteration for inversion (short division)
			
			int digits = (a_bits >= b_bits) ? a_bits / 4 : b_bits / 4;
			
			Number invb = b .set_precision(16 + digits) .inverse();
			
			//  Compute the unsigned quotient
			
			q = a .abs() .multiply(invb) .to_integer();
			
			Number remainder = a .subtract(q.multiply(b));
			
			//  Correct for a one-off error
			
			if (remainder.compare(b) >= 0) q = q .add(1);
			
			method = 4;
		}
		
		//  Replace the sign
		
		if (sign == '-') q = q .negate();
		
		
		
		//  Verify that the integer quotient is correct
		
		Number qxv = q.multiply(divisor);
		
		Number remainder = this->subtract(qxv);
		
		
		//  Compare the magnitudes of the remainder and divisor
		
		if ( ( remainder.abs() .compare(divisor.abs()) >= 0 ) ||
		
		     ( (this->sign != '-') && (remainder.sign == '-') ) ||
		     ( (this->sign == '-') && (remainder.sign != '-') && !remainder.equals(0) ) )
		{
			cout << "dividend == " <<   this->to_string(16) << endl;
			cout << "divisor  == " << divisor.to_string(16) << endl;
			cout << "quotient == " <<       q.to_string(16) << endl;
			cout << "this - q x v == " << remainder.to_string(16) << endl;
			cout << "division method == " << method << endl;
			
			string message = "\ninteger division error";
			
			cout << message << endl; throw string(message);
		}
		
		q = q .trim();
		
		return set_result(q);
	}
	
	
	//  Floating point division
	
	
	a = a.trim();  b = b.trim();
	
	int a_length, b_length;
	
	a_length = a .vec_int.size();
	b_length = b .vec_int.size();
	
	
	//  Set a and b equal to the same precision
	
	a = a .set_precision(precision);
	b = b .set_precision(precision);
	
	//  Calculate the difference between the two intpoints
	
	int a_b_intpoint = a.intpoint - b.intpoint;
	
	
	//  The dividend has to be left shifted so that the size of the
	//  quotient (== dividend length - divisor length) is greater than
	//  the precision of the dividend or divisor (whichever is greater).
	//  
	//  If the dividend length >= the divisor length, the dividend has to
	//  be left shifted by 4 bits x the precision (which should be rounded
	//  to a multiple of 32). The intpoint of the quotient is then set equal
	//  to the number of ints that were appended to the dividend.
	//
	//  If the dividend length < the divisor length, then it first has
	//  to be left shifted by an additional (dividend bits - divisor bits)
	//  to equalize the two lengths. This multiplies the quotient by a
	//  power of two so the quotient has to be right shifted after
	//  division to reduce it to the correct size. Right shifting is done
	//  by incrementing the intpoint of the quotient.
	//
	//  (To simplify the method, the dividend is left shifted by a multiple
	//  of 32 (or the dividend - divisor bits + 31 / 32) so that the integer
	//  point can be incremented instead of right shifting.)
	
	
	
	int left_shifts = 0;
	
	if (a.vec_int.size() < b.vec_int.size())
	{
		//  Equalize the two lengths
		
		int ints = b.vec_int.size() - a.vec_int.size();
		
		a = a .shift_left(32*ints, 32*ints);
		
		left_shifts += ints;
	}
	
	
	//  Left shift the dividend to the precision of the quotient
	
	int ints = 1 + precision / 8;  // 8 hex chars == 1 int
	
	a = a .shift_left(32*ints, 32*ints);
	
	left_shifts += ints;
	
	
	//  Compute the integer quotient and convert to floating point
	
	
	Number q;
	
	
	//  If the size of the divisor is much smaller than the dividend,
	//  the division operation can be done in parts by partitioning
	//  the dividend into smaller digits, computing the reduced
	//  quotients and remainders with the formulas
	//
	//  q[i+1] == ( r[i+0] + digit[i+0] ) * b / divisor, and
	//  r[i+1] == ( r[i+0] + digit[i+0] ) * b % divisor,
	//
	//  and then concatenating the quotients.
	
	
	
	//  Delete the intpoints of a and b to use an integer divider
	
	a.intpoint = 0;  a.precision = 0;
	b.intpoint = 0;  b.precision = 0;
	
	a.sign = '+';    b.sign = '+';
	
	
	
	if (false) {  }
	
	
//	else if ((a.vec_int.size() > 1) && (b.vec_int.size() == 1))
//	{
//		//  This method is not finished
//		
//		q = a.divide_by_int(b.int_value());
//		
//		//  Round the quotient up or down to the nearest integer
//		
//		if (q.greater_than(0))   q = q     .add(0.1);
//		else                     q = q.subtract(0.1);
//		
//		q = q .to_integer();
//	}
	
	
	else if ((a.vec_int.size() < 16) && (b.vec_int.size() < 16))
	{
		//  quadratic division (long division)
		
		q = a.quad_divide(b);
	}
	
	
	else
	{	//  Use the inverter to compute the quotient
		//
		//  a / b  ==  a x inv(b)
		
		//  First compute the inverse of the divisor
		
		Number inv_b = b.set_precision((a.vec_int.size() + 7) * 8).inverse();
		
		//  Multiply the inverse by the dividend
		
		q = a.multiply(inv_b);   q = q .to_integer();
	}
	
	
	
	//  Reduce the quotient by adding the left shifts to the intpoint and
	//  adding the difference between the dividend and divisor intpoints
	
	q.intpoint += (left_shifts + a_b_intpoint);
	
	//  If the signs are the same, then the quotient is positive
	//  If the signs are opposite, then the quotient is negative
	
	if (sign == '-') q = q .negate();
	
	q = q.set_precision(precision);
	
	
	//  Verify that the quotient is correct
	
	//  Compute the remainder of this - quotient x divisor
	
	Number qxv = q.multiply(divisor);
	
	Number remainder = this->subtract(qxv);
	
	if (remainder.abs() .compare(divisor.abs()) >= 0)
	{
		cout << "remainder == " << remainder.to_string(16) << endl;
		cout << "  divisor == " <<   divisor.to_string(16) << endl;
		
		cout << endl;
		
		cout << "remainder is greater than divisor" << endl;
		
		string message = "\nfloating point division error";
		
		cout << message << endl; throw string(message);
	}
	
	q = q.set_precision(precision);
	
	return set_result(q);
}



Number& Number::divide_by_int(int divisor)
{
	//  divides a large number by an int in O(n) operations
	
	//  ...
	
	//  ...
}



double Number::double_value()
{
	//  returns the value of this as a double
	
	Number number (*this);
	
	if (number.equals(0)) return 0;
	
	int exp = 0;
	
	//  sign bit of vec_int[0] should be zero
	
	if ((number.vec_int[0] & 0x80000000) != 0)
	{
		number = number.shift_right(1);
		
		exp += 1;
	}
	
	else
	{	while ((number.vec_int[0] & 0x40000000) == 0)
		{
			number = number.shift_left(1);
			
			exp -= 1;
		}
	}
	
	int int1 = number.vec_int[0], int2 = 0;
	
	if (number.vec_int.size() > 1)
	
	    int2 = number.vec_int[1];
	
	exp += (number.vec_int.size() - 1) * 32;
	
	exp -= number.intpoint * 32;
	
	double d = Math::pow(2, exp) *
	
	    (int1 + int2 * 1.0 / 0x10000 / 0x10000);
	
	if (this->less_than(0)) d *= -1;
	
	return d;
}


Number& Number::e_()
{
	return e_(32);
}

Number& Number::e_(int digits)
{
	int k = 128 + 256 * digits / 4096;
	
	int precision = digits + 4;
	
	Number x = Number(2).pow(k)
	
	    .set_precision(precision);
	
	x = x.inverse();
	
	Number e(0);
	
	//  iterations == digits / digits per iteration
	
	int t = 4 + digits / (k / 4);
	
	Number temp(1);
	
	for (int i = 1; i < t; i++)
	{
		e = e .add(temp);
		
		temp = temp .multiply(x).divide(i);
		
		//  cout << "e^(1/2^k) == " << e << endl;
	}
	
	e = e .pow(Number(2).pow(k));
	
	e = e .set_precision(digits);
	
	return *new Number(e);
}



bool Number::equals(int val)
{
	Number a = (this->trim());
	
	if (a.vec_int.size() > 1) return false;
	
	if (a.vec_int.size() == 1)
	{
		if (a.sign != '-')
		
		    return a.vec_int[0] == + val ? true : false;
		
		else if (a.sign == '-')
		
		    return a.vec_int[0] == - val ? true : false;
	}
	
	
	Number b (val);
	
	vector<int> vec_int1 = a.vec_int;
	vector<int> vec_int2 = b.vec_int;
	
	//  The number zero is allowed to be plus or minus
	//
	//  because multiplying zero by minus one is still zero
	
	if ((a.sign != b.sign) && (val != 0)) return false;
	
	return ( Math::compare(vec_int1, vec_int2) == 0 ) ? true : false ;
}


bool Number::equals(long val)
{
	Number number(val);
	
	return this->equals(number);
}

bool Number::equals(double val)
{
	//  Calculate the difference and compare to 0.000...0001
	
	Number a = (this->trim());
	
	Number b = a.subtract(val);
	
	return (b .abs() .double_value() < 0.00000001) ? true : false;
}

bool Number::equals(Number& val)
{
	//  tests for equality using the
	//  smaller of the two precisions
	
	return this->compare(val) == 0;
}


Number& Number::exp(int x)
{
	Number number(x);
	
	return exp(number);
}

Number& Number::exp(Number& x)
{
	//  returns e ^ x
	
	int precision = 16;
	
	if (x.precision > precision)
	
	    precision = x.precision;
	
	return e_(precision) .pow(x);
}


Number& Number::flip_bit(long bit)
{
	//  flips a bit
	
	Math::flip_bit(this->vec_int, bit);
	
	return *this;
}


float Number::float_value()
{
	//  returns a float
	
	double d = this->double_value();
	
	return (float) d;
}


Number& Number::gcd(long n)
{
	Number number(n);
	
	return this->gcd(number);
}

Number& Number::gcd(Number& n)
{
	//  returns the greatest common divisor
	
	//  The gcd function is only defined for integers
	
	if (!this->is_integer() || !n.is_integer())
	{
		string message = "gcd is only defined for integers";
		
		cout << message << endl; throw string(message);
	}
	
	Number a, b, c, m (*this);
	
	if (m.abs().greater_than(n.abs()))
	{
		a = m.abs();
		b = n.abs();
	}
	
	else
	{	a = n.abs();
		b = m.abs();
	}
	
	if (a.equals(0) || a.equals(1)
	 || b.equals(0) || b.equals(1))
	{
		return set_result(Number(1).lvalue());
	}
	
	if (a.equals(b))
	{
		return set_result(a);
	}
	
	//  Compute the gcd
	
	Number gcd(0);
	
	while (true)
	{
		c = a.mod(b);
		
		if (c.equals(0))
		
		    return (!gcd.equals(0)) ?
		
			set_result(gcd) : b;
		
		else { gcd = c; a = b; b = c; }
	}
}


int Number::get_bit(long bit)
{
	//  returns a bit
	
	int length = this->vec_int.size();
	
	if ((1 + bit / 32) > length) return 0;
	
	return Math::get_bit(this->vec_int, bit);
}


long Number::get_lowest_set_bit()
{
	return Math::get_lowest_set_bit(this->vec_int);
}

int Number::get_precision()
{
	//  returns the precision
	
	return this->precision;
}


int Number::int_value()
{
	//  returns the lowest 32-bit int
	
	Number n (*this);
	
	if (n.intpoint != 0) n = n.to_integer();
	
	int vector_length = n.vec_int.size();
	
	int intval = n.vec_int[vector_length - 1];
	
	if (n.sign == '-') intval = ~intval + 1;
	
	return intval;
}


Number& Number::inverse()
{

	//  uses the quadratically convergent iterative formula
	//  u = u (2 - v u) to compute the inverse u == 1 / v
	//
	//  This iteration requires only a few multiprecision multiplications
	
	
	//  The precision expands by 2 x the number of integer digits
	//
	//  The precision = 2 x int digits + fraction digits.
	//
	//  For example, the output of
	//
	//  new Number("12345678.0000").inverse());
	//
	//  should be  0.00000008100000664200
	//
	//  (because 2 x 8 + 4 == 20 digits of precision)
	
	
	if (this->is_complex())
	{
	
		//  Compute the reciprocal of a complex number a + i b
		//
		//         1     (a - i b)       a - i b
		//  ==  ________  ________  ==  _________
		//
		//      a + i b  (a - i b)      a^2 + b^2
		
		
		Number d = this->to_real().square()
		      .add(this->to_imag().square());
		
		Number real = this->to_real()         .divide(d);
		Number imag = this->to_imag().negate().divide(d);
		
		Number complex (real, imag);
		
		return set_result(complex);
	}
	
	
	//  Define zero (or the precision of zero)
	
	Number zero ("0", 16);
	
	zero = zero .set_precision(this->precision);
	
	if (this->equals(zero))
	{
		string message = "divide by zero";
		
		throw string(message);
	}
	
	if (this->bit_count() < 32)
	{
		Number inv (1.0 / this->vec_int[0]);
		
		return set_result(inv);
	}
	
	//  Calculate the number of int and frac digits
	
	int  intdigits = this->to_integer() .to_string(16).length();
	int fracdigits = this->to_fraction().to_string(16).length();
	
	int inv_precision = 2 * intdigits + fracdigits + 8;
	
	
	//  Compute u = 1 / v by iterating
	//
	//  u = u (2 - v u);  log(n) times
	
	
	//  Set u approximately equal to 1 / v.
	//
	//  Make sure the product v u is close to one
	//  or else u will converge to the inverse of v
	//  in O(n) instead of O(log(n)) multiplications
	
	
	Number v = (this->trim()), v1;
	
	if (v.precision == 0)  v = v.set_precision(8);
	
	
	//  Move the integer point of the input so that v ~ 1
	
	//  integer bits == v bits - 32 * intpoint
	
	const int d = (int) (v.bit_count() - 32 * v.intpoint);
	
	const int bits = (d / 32) * 32;
	
	if (d >= 0)
	{
		v1 = v .shift_left(bits + 32, bits);
		
		v1.intpoint += bits / 32;
		
		v1 = v1 .shift_right(d).trim();
	}
	
	else // if (d < 0)
	{
		v1 = v .shift_left(-d + 32, -d);
	}
	
	
	//  Set the initial value of u1 ~ 1 / v1
	
	const double dblv1 = v1.set_precision(16) .double_value();
	
	const double dblu1 = 1.0D / dblv1;
	
	Number u, u1 ((long) dblu1);
	
	
	
	//  First iterate at 32-bits precision
	
	int p = 8;
	
	for (int i = 0; i < 8; i++)  // constant precision
	{
		Number v1u1 = v1.set_precision(8) .multiply(u1.set_precision(p));
		
		Number two = Number(2) .set_precision(8);
		
		u1 = u1.set_precision(p) .multiply(two.subtract(v1u1));
	}
	
	
	//  Do a few more iterations doubling the precision for each iteration
	
	for (  ; (p < 1 * inv_precision) || (p < 64);  p *= 2)
	{
		Number v1u1 = v1.set_precision(p) .multiply(u1.set_precision(p));
		
		Number two = Number(2) .set_precision(p);
		
		u1 = u1.set_precision(p) .multiply(two.subtract(v1u1));
	}
	
	
	//  Move the integer point of the inverse
	
	if (d >= 0)
	{
		u1 = u1 .shift_left(bits, bits);
		
		u1 .intpoint += bits / 32;
		
		u = u1 .shift_right(d).trim();
	}
	
	else    u = u1 .shift_left(-d + 32, -d);
	
	
	//  Verify the inverse
	
	//  ...    ...
	
	
	u = u .set_precision(inv_precision) .trim();
	
	return set_result(u);
}


bool Number::is_base_64(string &str)
{
	//  tests if a string is encoded in base 64
	
	bool bool1 = true;
	
	char *array1 = Convert::string_to_char_array(str);
	
	int array_length1 = str.length();
	
	for (int i = 0; i < array_length1; i++)
	{
		const int *intarray = Convert::base_64_to_int;
		
		int array_length = sizeof(Convert::base_64_to_int)
		                 / sizeof(Convert::base_64_to_int[0]);
		
		if (array1[i] >= array_length)
		
		    { bool1 = false; break; }
		
		if ((Convert::base_64_to_int[
		
		    array1[i]] == -1) && (array1[i] != '='))
		
		    { bool1 = false; break; }
	}
	
	if ((str.length() % 4) != 0)
	
	    bool1 = false;
	
	return bool1;
}


bool Number::is_complex()
{
	//  tests if the number is complex (a + b i)
	
	if (this->vec_int1.empty()) return false;
	
	Number number(this->vec_int1);
	
	Number zero = Number(0) .set_precision(this->precision1);
	
	if (number.equals(zero)) return false;
	
	return true;
}

bool Number::is_coprime_with(long n)
{
	Number number(n);
	
	return this->is_coprime_with(number);
}

bool Number::is_coprime_with(Number& n)
{
	Number number(*this);
	
	return number .gcd(n) .equals(1) ? true : false;
}

bool Number::is_digit_string(string &str, int radix)
{
	if (str.find("+") != string::npos || str.find("-") != string::npos)
	
	    return false;
	
	return is_integer_string(str, radix);
}

bool Number::is_divisible_by(int n)
{
	return this->mod(n).equals(0);
}

bool Number::is_divisible_by(Number& n)
{
	return this->mod(n).equals(0);
}

bool Number::is_even()
{
	return this->get_bit(0) == 0 ? true : false;
}


bool Number::is_generator(Number& p)
{
	//  tests whether a number is a generator modulo p
	
	if (!p.is_prime()) throw string("non-prime modulus");
	
	int limit = 1024 * 1024;
	
	vector<int> factors = Math::factor(p.subtract(1), limit);
	
	int factor_length = factors.size();
	
	for (int i = 0; i < factor_length; i++)
	{
		int factor = factors[i];
		
		Number residue = this ->mod_pow(
		
		    p.subtract(1) .divide(factor), p);
		
		if (residue .equals(1)) return false;
	}
	
	return true;
}


bool Number::greater_than(double val)
{
	Number number(val);
	
	return this->greater_than(number);
}

bool Number::greater_than(int val)
{
	Number number(val);
	
	return this->greater_than(number);
}

bool Number::greater_than(long val)
{
	Number number(val);
	
	return this->greater_than(number);
}

bool Number::greater_than(Number& val)
{
	if (this->equals(val)) return false;
	
	else return this->compare(val) == 1;
}

bool Number::is_integer()
{
	// returns true if the precision is zero
	
	Number n (*this);
	
	return (n.precision  == 0)
	    && (n.precision1 == 0) ? true : false;
}

bool Number::is_integer_string(string& str, int radix)
{
	if (str.empty()) return false;
	
	if ( (str.find(" ") != string::npos)
	  || (str.find(".") != string::npos) )
	
	    return false;
	
	return is_number_string(str, radix);
}

bool Number::less_than(double val)
{
	Number number(val);
	
	return this->less_than(number);
}

bool Number::less_than(int val)
{
	Number number(val);
	
	return this->less_than(number);
}

bool Number::less_than(long val)
{
	Number number(val);
	
	return this->less_than(number);
}

bool Number::less_than(Number& val)
{
	Number n (*this);
	
	return n.compare(val) == -1;
}

bool Number::is_number_string(string& str, int radix)
{
	//  A number string can only have a string of characters
	//  in base radix, an intpoint, and a sign.
	
	//  An example of a number string, integer string, and
	//  digit string in base-16
	//
	//  - 1.23456789abcdef //  number string (chars, sign, and intpoint)
	//   - 123456789abcdef // integer string (chars and sign)
	//     123456789abcdef //   digit string (chars only)
	
	if (str.empty()) return false;
	
	string s = string(str);  int pos = 0;
	
	while ((pos = s.find(" ")) != string::npos)
	
	    s.erase(pos, 1);
	
	
	//  If string starts with a sign, remove the sign (and the space)
	
	if ((s.find_first_of("+") == 0) || (s.find_first_of("-") == 0))
	{
		s = s.substr(1, s.length() -1);
		
		if (s.find_first_of(" ") == 0)
		
		    s = s.substr(1, s.length() -1);
	}
	
	if (s.empty()) return false;
	
	//  only one integer / fraction point is allowed
	
	if (s.find_first_of('.') != s.find_last_of ('.'))
	
	    return false;
	
	//  erase the integer / fraction point
	
	if (s.find('.') != string::npos)
	
	    s = s .erase(s.find('.'), 1);
	
	for (int i = 0; i < s.length(); i++)
	{
		int digit = 0;
		
		char c = str[i];
		
		if  (isalpha(c)) digit = c - 'a' + 10;
		else             digit = c - '0' +  0;
		
		//  cout << "char digit == " << c << "  " << digit << endl;
		
		if (digit >= radix) return false;
	}
	
	return true;
}


bool Number::is_power(int n)
{
	Number number (*this);
	
	if (number.length() <= 2)
	
	    return Math::is_power(number.long_value(), n);
	
	Number root = number.root(n).add(0.1).to_integer();
	
	return root.pow(n) .equals(number) ? true : false;
}


bool Number::is_power_of(int n)
{
	Number number (*this);
	
	if (number.equals(1) || number.equals(n)
	
	 || number.equals(Number(n).square()))
	
	    return true;
	
	
	//  for base == 2
	
	//  2^k & (2^k -1) == 10000... & 01111... == 0
	
	if (n == 2)  return number.and_(
	
	    number.subtract(1)).equals(0) ? true : false;
	
	
	//  Compute number (mod base)
	
	if (!number.is_integer() || number.equals(0)
	
	 || !number.mod(n).equals(0))
	
	    return false;
	
	
	//  Reduce the size of the dividend
	
	long a_size = number.bit_count();
	
	long b_size = Math::bit_count(n);
	
	int d_size = (int) (a_size - b_size) / Math::log2(n);
	
	Number divisor = Number(n) .pow(d_size);
	
	Number quotient = number.divide(divisor);
	
	if (quotient.equals(1)) return true;
	
	while (quotient.greater_than(1))
	{
		quotient = quotient.divide(n);
		
		if (quotient.equals(1))
		
		    return true;
	}
	
	return false;
}


bool Number::is_prime()
{
	return is_probable_prime(2);
}


bool Number::is_probable_prime()
{
	return this->is_probable_prime(2);
}

bool Number::is_probable_prime(int certainty)
{

	Number number (*this);
	
	if (!number.is_integer())
	{
		string message = "number must be an integer";
		
		throw string(message);
	}
	
	Number n = number.trim();
	
	
	if (n.equals(0) || n.equals(1)) return false;
	
	if (n.equals(2) || n.equals(3)) return true;
	
	if (n.is_even() && !n.equals(2)) return false;
	
	int primes[168]
	{
	     2,    3,    5,    7,   11,   13,   17,   19,   23,   29,   31,   37,
	    41,   43,   47,   53,   59,   61,   67,   71,   73,   79,   83,   89,
	    97,  101,  103,  107,  109,  113,  127,  131,  137,  139,  149,  151,
	   157,  163,  167,  173,  179,  181,  191,  193,  197,  199,  211,  223,
	   227,  229,  233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
	   283,  293,  307,  311,  313,  317,  331,  337,  347,  349,  353,  359,
	   367,  373,  379,  383,  389,  397,  401,  409,  419,  421,  431,  433,
	   439,  443,  449,  457,  461,  463,  467,  479,  487,  491,  499,  503,
	   509,  521,  523,  541,  547,  557,  563,  569,  571,  577,  587,  593,
	   599,  601,  607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
	   661,  673,  677,  683,  691,  701,  709,  719,  727,  733,  739,  743,
	   751,  757,  761,  769,  773,  787,  797,  809,  811,  821,  823,  827,
	   829,  839,  853,  857,  859,  863,  877,  881,  883,  887,  907,  911,
	   919,  929,  937,  941,  947,  953,  967,  971,  977,  983,  991,  997,
	};
	
	int primes_length = sizeof(primes) / sizeof(primes[0]);
	
	int bits = (int) n.bit_count();
	
	if (number.bit_count() < 16)
	
	    for (int i = 0; i < primes_length; i++)
	
		if (number.equals(primes[i])) return true;
	
	//  Divide by small primes
	
	for (int i = 0; i < primes_length; i++)
	
	    if (number.is_divisible_by(primes[i])) return false;
	
	//  Apply the Fermat test a^(n-1) mod n == 1
	//  using bases a = 2 and a = 3 as witnesses
	
	//  If n is prime, then the totient is phi(n) == n-1.
	
	if ( !Number(2) .mod_pow(n.subtract(1), n) .equals(1)
	  || !Number(3) .mod_pow(n.subtract(1), n) .equals(1) )
	
		return false;  // composite
	
	
	//  If the number has made it this far, it has passed the Fermat test
	//  using the first two primes as bases.
	//
	//  A composite number n that can pass the Fermat test a^(n-1) == 1 (mod n)
	//  for some base a is called a pseudo prime. A composite number that can
	//  pass the Fermat test for all bases is called an absolute pseudo prime
	//  or Carmichael number.
	
	
	//  Test for Carmichael numbers and other pseudo primes
	
	//  cout << "Testing for pseudo primes" << endl;
	
	
	//  n - 1 = (2 ^ s) r
	
	Number r = number.subtract(1);  int s = 0;
	
	while (r.is_even()) { r = r.shift_right(1); s++; }
	
	for (int i = 1; i <= certainty; i++)
	{
		//  Choose 1 < a < n-1
		
		Number a (i + 1),  y;
		
		//  Compute y = a^r (mod n) and test if it equals 1 or -1
		
		if ( !(y = a.mod_pow(r, n)) .equals(1) && !y .equals(n.subtract(1)) )
		{
			int j = 1;
			
			while ( (j++ <= (s - 1)) && !y.equals(n.subtract(1)) )
			{
				//  Compute y = y ^ 2  mod n;
				
				y = y .square() .mod(n);
				
				if (y.equals(1)) return false;
			}
			
			if (!y.equals(n.subtract(1))) return false;
		}
	}
	
	return true;
}



Number& Number::lcm(long n)
{
	Number number(n);
	
	return this->lcm(number);
}

Number& Number::lcm(Number& a)
{
	//  returns the least common multiple
	
	Number number (*this);
	
	Number lcm = number.multiply(a)
	
	    .divide(number.gcd(a));
	
	return set_result(lcm);
}


int Number::length()
{
	//  returns the length in integers
	
	Number number (*this);
	
	return number.vec_int.size();
}


Number& Number::ln()
{
	//  computes the natural logarithm
	//
	//  for any real number x > 0
	
	Number number (*this);
	
	Number ln = number.log();
	
	return set_result(ln);
}

long Number::log2()
{
	//  returns the bit count
	
	//  To compute the log to the base 2
	//  for a real number n, use n.log(2)
	
	Number number (*this);
	
	return number.to_integer() .bit_count();
}

Number& Number::log(int base)
{
	//  log(n) == log(n) / log(b)
	//     b         x        x
	
	Number Base (base);
	
	Number number (*this);
	
	Number log = number.log(Base);
	
	return set_result(log);
}

Number& Number::log(Number& base)
{
	//   To compute the log to any base b, use the formula
	//
	//   log n  ==  log n / log b
	//      b          x       x
	
	Number number (*this);
	
	Number log = number.log().divide(base.log());
	
	return set_result(log);
}

Number& Number::log()
{
	//  computes the natual log of a number
	
	Number number (*this);
	
	if (number.equals(0))
	{
		string message = "log(0) == - infinity";
		
		cout << message << endl; throw string(message);
	}
	
	if (number.less_than(0))
	{
		//  from  e^(i pi) == -1, log(-1) == i pi
		
		string message = "log of -x == log(-1) + log(x) == i pi + log(x)";
		
		cout << message << endl; throw string(message);
	}
	
	Number arg(number);
	
	if (arg.precision < 16)
	
	    arg.precision = 16;
	
	
	//  Move the integer / fraction point so the argument is close to 1
	
	int exp = 0;
	
	while (arg.greater_than(0x7fffffff)) { arg.intpoint++; exp += 32; }
	
	while (arg.greater_than(1)) { arg = arg.shift_right(1); exp++; }
	
	while (arg.intpoint > arg.vec_int.size()) { arg.intpoint--; exp -= 32; }
	
	while (arg.less_than(1)) { arg = arg.multiply(2); exp--; }
	
	if (arg.greater_than(1.8)) { arg = arg.divide(2); exp++; }
	
	
	const int p = 2 * arg.get_precision();
	
	
	//  Compute log(arg x 2^exp) == log(arg) + exp log(2)
	
	//  The argument of log(2) is too large so first we divide by e
	//
	//  log(2 / e) == log(2) - log(e) == log(2) - 1
	//
	//  Therefore log(2) == log(2 / e) + 1
	
	Number arg1 = arg.set_precision(precision);
	
	Number log1 = arg1 .log1();
	
	Number log2 = Number(2) .divide(
	
	    e_(precision)) .log1() .add(1);
	
	Number log = log1 .add(log2.multiply(exp));
	
	
	//  Check the answer
	//
	//  e ^ log(x) / x == 1
	
	double e2logx = Math::pow(Math::e, log.double_value());
	
	//  System.out.println(e2logx / this->double_value());
	
	
	return set_result(log);
}


Number& Number::log1()
{

	//  The logarithmic function log(1 + x) is computed from the
	//
	//             n+1   n
	//  series (-1)    x  / n  which is defined
	//
	//  for -1 < x <= 1  or  0 < 1 + x <= 2.
	//
	//  (The argument 1 + x is used instead of x
	//
	//  because log(x = 0) is undefined)
	//
	//  The argument (1 + x) has to be close to 1
	//
	//  or the series will not converge
	
	if (this->greater_than(2) || this->less_than(0))
	{
		string message = "argument is not between 0 and 2";
		
		throw string(message);
	}
	
	Number log(0);
	
	Number number (*this);
	
	Number x = this->subtract(1);
	
	int bits = 4 + this->precision;
	
	Number powerofx = Number(1) .set_precision(bits);
	
	for (int n = 1;  ; n++)
	{
		//                    n+1   n
		//  log(1 + x) == (-1)    x  / n  for 0 < (1 + x) < 2
		//
		//  log = log.add( x.pow(n) .divide(n)
		//
		//  	.multiply(((n % 2) == 1) ? 1 : -1));
		
		powerofx = powerofx .multiply(x);
		
		Number log1 = powerofx .divide(n) .multiply(((n % 2) == 1) ? 1 : -1 );
		
		log = log .add(log1);
		
		if (log1.equals(0))  break;
	}
	
	return set_result(log);
}


long Number::long_value()
{
	Number n (*this);
	
	if (n.intpoint != 0) n = n.to_integer();
	
	int nvector_length = n.vec_int.size();
	
	long longval = n.vec_int[nvector_length - 1];
	
	longval &= 0xffffffffL;
	
	if (nvector_length > 1)
	
	    longval += ((1L * n.vec_int[nvector_length - 2]) << 32);
	
	longval &= 0x7fffffffffffffffL;
	
	if (n.sign == '-')
	
	    longval = ~longval + 1;
	
	return longval;
}


Number& Number::lvalue()
{

	//  This method converts an rvalue to an lvalue so that it can
	//  be used to initialize the constructor of a Number because
	//  the compiler will complain if the user tries to initialize
	//  a Number constructor using an rvalue.
	//
	//  The difference between lvalues and rvalues is that lvalues
	//  are persistent and can be assigned to other values. rvalues
	//  are ephemeral and cannot be assigned to other values. (In
	//  the C language lvalue and rvalue mean left and right value.)
	//
	//  For instance, Number(0) is an rvalue because it represents
	//  the value of zero but has no variable name and cannot be
	//  assigned to another value. But Number number(0) is an lvalue
	//  because it can be used on the left side of an equation and it
	//  can be assigned to another value.
	//
	//  The lvalue method doesn't return a variable name but it does
	//  return an anonymous Number object created on the heap using
	//  the new keyword which makes it an lvalue. This object will
	//  persist until the compiler deletes it by calling the destructor
	//  when the object goes out of scope.
	//
	//  The following code compiles on the g++ compiler and doesn't
	//  cause a memory leak which proves that the compiler is deleting
	//  the objects created by the lvalue method.
	//
	//  for (int i = 0; i < 1234567890; i++)
	//
	//      Number number(Number(0).lvalue());
	
	
	return set_result(*this);
}


Number& Number::max(Number& val)
{
	//  returns the greater of two numbers
	
	return (this->compare(val) >= 0) ? *this : val;
}

Number& Number::min(Number& val)
{
	//  returns the lesser of two numbers
	
	return (this->compare(val) < 0) ? *this : val;
}

Number& Number::mod(double val)
{
	Number number(val);
	
	return this->mod(number);
}

Number& Number::mod(int val)
{
	//  This method is faster than the mod(number) method
	//  because it runs in O(n) time instead of O(n^2)
	
	//  Compute the mods of the ints and add them
	
	long residue = 0L;
	
	long twosr = 1L; // twos residue
	
	int this_length = this->length();
	
	for (int i = this_length - 1; i >= 0; i--)
	{
		long a = this->vec_int[i] & 0xffffffffL;
		
		long r = a % val;
		
		r = (r * twosr) % val;
		
		twosr *= 0x100000000L;
		
		twosr %= val;
		
		residue += r;
	}
	
	residue %= val;
	
	return set_result(Number(residue).lvalue());
}

Number& Number::mod(long val)
{
	Number number(val);
	
	return mod(number);
}

Number& Number::mod(Number& mod)
{
	//  reduces a number a modulo n
	
	//  Note that this method works similarly to the int operator
	//  % which reduces the number to the range - n < residue < + n.
	//
	//  If the argument a is negative but the residue r has to be positive
	//  then use r = a .mod(n) .add(n) .mod(n) or r = ((a % n) + n) % n.
	
	
	if (mod.is_complex())
	{
		string message = "complex modulus";
		
		cout << message << endl; throw string(message);
	}
	
	if (this->is_complex())
	{
		Number real = this->to_real();
		Number imag = this->to_imag();
		
		real = real .mod(mod) .add(mod) .mod(mod);
		imag = imag .mod(mod) .add(mod) .mod(mod);
		
		Number complex (real, imag);
		
		return set_result(complex);
	}
	
	if (mod.equals(0))
	{
		string message = "divide by zero";
		
		cout << message << endl; throw string(message);
	}
	
	//  The modulus cannot be negative or the program will hang
	
	Number a = Number(*this).abs();
	Number n = Number(mod)  .abs();
	
	char sign = (this->sign == '+') ? '+' : '-';
	
	if (a.less_than(n))
	{
		a.sign = sign;
		
		return set_result(a);
	}
	
	if (a.equals(0))
	
	    return set_result(Number(0).lvalue());
	
	//  Compute the residue
	//
	//  r = a (mod n) = a - [a / n] n
	
	Number q = a.divide(n);
	
	Number qxn = q .to_integer() .multiply(n);
	
	Number r = a.subtract(qxn);
	
	r.sign = sign;
	
	return set_result(r);
}


Number& Number::mod(Number& n, Number& inv)
{


}


Number& Number::mod_divide(int divisor, Number& n)
{
	Number number (divisor);
	
	return this->mod_divide(number, n);
}

Number& Number::mod_divide(Number& divisor, Number& n)
{

	//  modular division
	
	//  This method allows the user to compute a / b  (mod n)
	//
	//  using  a .modDivide(b, n)  instead of
	//
	//  a .multiply( b.modInverse(n) ) .mod(n)
	
	
	if (!this->is_integer() || !divisor.is_integer() || !n.is_integer())
	{
		string message = "modular division is only defined for integers";
		
		cout << message << endl; throw string(message);
	}
	
	Number a = this->mod(n);
	
	if (divisor.equals(0) || !divisor.is_coprime_with(n))
	{
		string message = "divisor and modulus are not coprime";
		
		cout << message << endl; throw string(message);
	}
	
	//  Compute the modular quotient
	
	Number quotient = a .multiply(divisor.mod_inverse(n)) .mod(n);
	
	if (!quotient .multiply(divisor) .mod(n) .equals(a))
	{
		string message = "mod divide error";
		
		cout << message << endl; throw string(message);
	}
	
	return set_result(quotient);
}


Number& Number::mod_inverse(long modulus)
{
	Number number(modulus);
	
	return this->mod_inverse(number);
}

Number& Number::mod_inverse(Number& n)
{
	//  computes the modular inverse of a number
	
	if (this->is_complex())
	{
		Number d = this->to_real().square().mod(n)
		      .add(this->to_imag().square().mod(n)) .mod(n);
		
		Number real = this->to_real()          .mod_divide(d, n);
		Number imag = this->to_imag().negate(n).mod_divide(d, n);
		
		Number complex (real, imag);
		
		return set_result(complex);
	}
	
	if (!this->is_integer() || !n.is_integer())
	{
		string message = "modular inversion is only defined for integers";
		
		cout << message << endl; throw string(message);
	}
	
	if (this->equals(1))
	
		return set_result(*this);
	
	char sign = (this->sign == '+') ? '+' : '-';
	
	if (this->equals(0) || !this->is_coprime_with(n))
	{
		string message = "non-invertible number";
		
		cout << message << endl; throw string(message);
	}
	
	//  Use the extended Euclidean algorithm
	
	Number a = Number(*this).abs();
	
	Number inva = a .mod_inverse1(n);
	
	
	//  Replace the sign
	
	if (sign == '-')
	
	    inva = n.subtract(inva);
	
	
	//  Check the answer
	
	if (!inva .multiply(*this).mod(n).add(n).mod(n).equals(1))
	{
		string message = "modular inversion or division error";
		
		cout << message << endl; throw string(message);
	}
	
	return set_result(inva);
}


Number& Number::mod_inverse1(Number& n)
{
	//  computes the modular inverse of a number
	
	//  This method implements the extended Euclidean algorithm
	
	//  variables:  a, b, c, q, r,   x, x1, x2,   y, y1, y2;
	
	Number a, b, c, q, r;
	
	Number x, x1, x2;
	Number y, y1, y2;
	
	x1 = Number(0).lvalue();  x2 = Number(1).lvalue();
	y1 = Number(1).lvalue();  y2 = Number(0).lvalue();
	
	a = this ->abs() .mod(n);
	
	b = n .abs();
	
	while (!a.equals(0))
	{
		q = b.divide(a);
		
		r = b.subtract(q.multiply(a));
		
		x = x2.subtract(q.multiply(x1));
		y = y2.subtract(q.multiply(y1));
		
		x2 = x1;  x1 = x;
		y2 = y1;  y1 = y;
		
		b = a;  a = r;
	}
	
	c = b; x = x2; y = y2;
	
	Number inva = y;
	
	if (!c.equals(1))
	{
		string message = string("number ") + this->to_string() +
		
		  string(" is not invertible ") + string("modulo ") + n.to_string(10);
		
		message += string("\n\ngcd =  ") + c.to_string();
		
		cout << message << endl; throw string(message);
	}
	
	inva.sign = sign;
	
	Number product = inva.multiply(*this).mod(n);
	
	if (product.abs() .equals(n.subtract(1)))
	
	    inva = inva.negate(n);
	
	return set_result(inva);
}


Number& Number::mod_multiply(Number& ar, Number& n, Number& n1, Number& r)
{

	//  Montgomery multiplication
	//
	//  This method performs modular multiplication without division.
	//  This is useful for methods that require large numbers of
	//  modular multiplications such as modular exponentiation.
	
	
	//  Montgomery multiplication
	//
	//  This method computes the product
	//
	//  zr == xr yr r^-1 (mod n)
	//
	//     == x y r^2 r^-1 == z r (mod n)
	//
	//  where
	//
	//  r == 2 ^ k;
	//
	//  n1 == - n ^-1 (mod r)
	//   == r - n ^-1 (mod r);
	//
	//  and (n, r) == 1.
	//
	//  r is chosen to be 2^k so that shifts and ands can be used
	//  instead of divs and mods (quotients and remainders).
	//
	//  In modular exponentiation, the r can be removed from the
	//  product zr by performing a final Montgomery multiplication
	//  with xr equal to zr, and yr equal to 1 (multiplication by 1);
	//  then (xr = zr) * (yr = 1) == zr * 1 * r^-1 == z.
	
	
	//  Compute  zr = xr * yr * r^-1 (mod n)
	
	
	Number xr = *this, yr = ar;
	
	int k = 32 * n.vec_int.size();
	
	Number rminus1 = r.subtract(1);
	
	Number xy;
	
	
	//  Compute z = ( xy + ( xy n1 mod r ) n ) / r
	//
	//            ==  xy  r^-1  (mod n)
	
	xy = xr.multiply(yr);
	
	Number temp1 = xy .and_(rminus1)
	
	    .multiply(n1) .and_(rminus1);
	
	Number temp2 = n.multiply(temp1);
	
	Number z = xy .add(temp2);
	
	if (k <= z.vec_int.size() * 32)
	
	     z = z.shift_right(k);
	
	else z = Number(0).lvalue();
	
	if (z.compare(n) >= 0)
	
	     z = z.subtract(n);
	
	return set_result(z);
}


Number& Number::mod_pow1(Number& exp, Number& n)
{

	//  Montgomery exponentiation
	
	//  Define
	//
	//  r = 2 ^ k,
	//
	//  n1 = - (n ^-1) mod r
	//
	//  a1 = a r mod(n)
	//  y1 = y r mod(n)
	//  x1 = x r mod(n)
	
	int k = 32 * n.vec_int.size();
	
	Number r = Number(2).pow(k);
	
	Number a = this->mod(n), y(1), x(exp);
	
	Number a1 = a .multiply(r) .mod(n);
	Number y1 = y .multiply(r) .mod(n);
	Number x1 = x .multiply(r) .mod(n);
	
	//  Compute  n' = r - n / r;
	
	Number n1 = r.subtract(n.mod_inverse(r));
	
	
	//  Compute y = a ^ x (mod n)
	
	//  Use the square and multiply method
	
	while (!x.equals(0))
	{
		//  Accumulate squares for each one-bit in the exponent
		
		if (x.test_bit(0))
		
		    y1 = a1 .mod_multiply(y1, n, n1, r);
		
		//  y' = a' y' = a r y r r^-1 mod n = a y r (mod n)
		
		//  Square the square
		
		a1 = a1 .mod_multiply(a1, n, n1, r);
		
		//  a' = a r a r r^-1 mod n = a^2 r (mod n)
		
		//  Shift the exponent to the next bit
		
		x = x.shift_right(1);
	}
	
	//  Remove the residue r from y1 = r a^x (mod n)
	
	Number one (1);
	
	y = y1 .mod_multiply(one, n, n1, r);
	
	return set_result(y);
}


Number& Number::mod_pow(int exp, int m)
{
	Number number1(exp), number2(m);
	
	return this->mod_pow(number1, number2);
}

Number& Number::mod_pow(int exp, Number& m)
{
	Number number1(exp), number2(m);
	
	return this->mod_pow(number1, number2);
}

Number& Number::mod_pow(Number& exp, int m)
{
	Number number1(exp), number2(m);
	
	return this->mod_pow(number1, number2);
}

Number& Number::mod_pow(Number& exp, Number& n)
{
	//  returns y = a ^ x (mod n)
	
	if (!exp.is_integer() || !n.is_integer())
	{
		string message = "";
		
		if (!exp.is_integer())  message =
		
		    "use the pow method for non-integer exponents";
		
		else if (!n.is_integer())  message =
		
		    "modular exponentiation is only defined for integers";
		
		cout << message << endl; throw string(message);
	}
	
	if (n.less_than(0))
	{
		string message = "negative modulus";
		
		cout << message << endl; throw string(message);
	}
	
	if (exp .equals(0))
	{
		return set_result(Number(1).lvalue());
	}
	
	if (this->equals(0))
	{
		return set_result(Number(0).lvalue());
	}
	
	
	//  Compute y = a^x (mod n)
	
	Number a = this->mod(n), y(1), x(exp);
	
	if (exp.less_than(0))
	{
		a = a.mod_inverse(n);
		
		x = x.abs();
	}
	
	//  if (!n.is_even())
	//  {
	//  	return a.mod_pow1(x, n);
	//  }
	
	
	//  Pre-compute the inverse for fast modular reduction
	
	const int digits = (int) n.bit_count() / 4;
	
	Number inv_n = n .set_precision(16 + digits) .inverse();
	
	
	//  Use the square and multiply method
	
	while (!x.equals(0))
	{
		//  Accumulate squares
		
		if (x.test_bit(0))
		{
			y = a.multiply(y);
			
			y = y.mod(n, inv_n);
		}
		
		//  Square the square
		
		a = a.square();
		
		a = a.mod(n, inv_n);
		
		//  Shift the exponent to the next bit
		
		x = x.shift_right(1);
	}
	
	y = y .mod(n);
	
	return set_result(y);
}


Number& Number::multiply(int multiplier)
{
	Number number(multiplier);
	
	return this->multiply(number);
}

Number& Number::multiply(long multiplier)
{
	Number number(multiplier);
	
	return this->multiply(number);
}

Number& Number::multiply(double multiplier)
{
	Number number(multiplier);
	
	return this->multiply(number);
}

Number& Number::multiply(Number& multiplier)
{
	//  multiplies two numbers using the Math multiply method
	
	//  Multiplication can also be done by squaring numbers using the identity
	//
	//  (a + b)^2 - (a - b)^2 == a^2 + 2 a b + b^2 - (a^2 - 2a b + b^2) == 4 a b
	//
	//  or  4 a b == the sum squared minus the difference squared.
	//
	//  (The 4 is removed by shifting two bits to the right)
	
	
	//  product == multiplicand x multiplier
	
	int precision = (this->precision > multiplier.precision) ?
	                 this->precision : multiplier.precision;
	
	Number product (Math::multiply(this->vec_int, multiplier.vec_int));
	
	product.intpoint = this->intpoint + multiplier.intpoint;
	
	product = product .trim() .set_precision(precision);
	
	product.sign = ((this->sign != multiplier.sign)
	
	    && !product.equals(0)) ? '-' : '+';
	
	//  This is used for debugging
	//
	//  double this_dbl    =      this->double_value();
	//  double mult_dbl    = multiplier.double_value();
	//  double product_dbl =    product.double_value();
	//  
	//  cout << "product this x multiplier == "
	//  
	//      << product_dbl << "  " << (this_dbl * mult_dbl) << endl;
	
	
	//  Return the product
	
	if (!this->is_complex() && !multiplier.is_complex())
	
	    return set_result(product);
	
	
	//  Compute the complex product
	//
	//     (a1 + i b1) * (a2 + i b2)
	//
	//  == (a1 a2 - b1 b2) + i (a1 b2 + a2 b1)
	
	
	Number real1, real2, imag1, imag2;
	
	real1 =       this->to_real();
	imag1 =       this->to_imag();
	
	real2 = multiplier.to_real();
	imag2 = multiplier.to_imag();
	
	Number real = real1.multiply(real2)
            .subtract(imag1.multiply(imag2));
	
	Number imag = real1.multiply(imag2)
	         .add(real2.multiply(imag1));
	
	Number complex (real, imag);
	
	return set_result(complex);
}


Number& Number::negate()
{
	Number n (*this);
	
	if      (this->sign == '+') n.sign = '-';
	else if (this->sign == '-') n.sign = '+';
	
	return set_result(n);
}

Number& Number::negate(int n)
{
	Number number(n);
	
	return this->negate(number);
}

Number& Number::negate(Number& n)
{
	//  negates the sign of a number
	
	Number number;
	
	if (this->greater_than(n))
	
	     number = n.subtract(this->mod(n));
	else number = n.subtract(*this);
	
	return set_result(number);
}


Number& Number::next_prime()
{
	// returns the next prime
	
	if (!this->is_integer())
	
	    throw string("illegal argument");
	
	Number n (*this);
	
	if (n.compare(2) == -1)
	
	    return set_result(Number(2).lvalue());
	
	if (n.is_even()) n = n.add(1);
	else             n = n.add(2);
	
	while (!n.is_prime())  n = n.add(2);
	
	return set_result(n);
}


int Number::parse_int(string &s)
{
	return parse_int(s, 10);
}


int Number::parse_int(string &s, int radix)
{
	if (!is_integer_string(s, radix)) return 0;
	
	Number number(s, radix);
	
	if (number.is_integer() && (number.compare((int) 2147483647) <= 0))
	
	     return number.int_value();
	
	else return 0;
}


Number& Number::pi()
{
	//  returns the value of pi
	
	return pi(32);
}

Number& Number::pi(int digits)
{
	//  computes the value of pi
	//  to any number of digits
	
	if (digits <= 256)
	
	      return pi11(digits);
	else  return pi1 (digits);
}


Number& Number::pi1(int digits)
{
	int d = digits + 64;
	
	Number a, y;
	
	a = Number(6) .set_precision(d) .subtract(
	    Number(4) .set_precision(d) .multiply(
	    Number(2) .set_precision(d) .sqrt() ) );
	
	y = Number(2) .set_precision(d) .sqrt() .subtract(1);
	
	for (int k = 0; k < Math::log2(d) / 2; k++)
	{
		//  k(d=16K) == Math::log2(16K)/2 ==  7 iterations
		//  k(d= 1M) == Math::log2(1M) /2 == 10 iterations
		
		//  or 1 quartic root + 1 cube + 5 squares + 1 divide
		//  + 2 multiplies == 16 multiplies per iteration
		
		Number z = Number(1) .set_precision(d) .subtract(y.pow(4)) .root(4);
		
		y = Number(1).set_precision(d) .subtract(z) .divide(
		    Number(1).set_precision(d) .add(z) );
		
		a = Number(1).set_precision(d) .add(y) .square().square()
		    .multiply(a) .subtract( y.cube() .add(y.square()) .add(y)
			.multiply(Number(2).pow(2*k + 3)) );
	}
	
	Number pi = a.inverse();
	
	pi = pi.set_precision(digits);
	
	return *new Number(pi);
}


Number& Number::pi11(int digits)
{

	//  This formula by Borwein and Borwein (B & B)
	//
	//  is twice as fast as Chudnovsky and Chudnovsky (C & C)
	//
	//                      n
	//   1        __    (-1) (6 n)! (A + B n)
	//  ____  ==  \     _____________________
	//            /_        3         n + 1/2
	//  12 pi     n=0   (n!) (3 n)! C
	//
	//                               __
	//  where     A = 212175710912 \/61 + 1657145277365
	//                                 __
	//            B = 13773980892672 \/61 + 107578229802750
	//                                         __   3
	//            C = [ 5280 (236674 + 30303 \/61) ]
	//
	//  Each term of this series adds 31 digits
	
	
	
	int d = digits * 6 / 5 + 64;
	
	Number c1 ("1657145277365");
	Number c2 ("107578229802750");
	
	Number A = Number("212175710912") .set_precision(d) .multiply(
	
	    Number(61).set_precision(d).sqrt()) .add(c1);
	
	Number B = Number("13773980892672") .set_precision(d) .multiply(
	
	    Number(61).set_precision(d).sqrt() ).add(c2);
	
	Number C = Number("5280") .set_precision(d) .multiply(
	
	    Number("236674") .add( Number("30303") .multiply(
	
		Number(61).set_precision(d).sqrt() ) ) ) .cube();
	
	Number D = C .sqrt() .inverse();
	
	Number D2 = D .square();
	
	Number E  = Number(1) .set_precision(d);
	Number F  = Number(1) .set_precision(d);
	Number G  = Number(1) .set_precision(d);
	Number H  = Number(1) .set_precision(d);
	Number H1 = Number(1) .set_precision(d);
	
	G = G .multiply(D);
	
	Number sum (0);
	
	Number product (1);
	
	for (int n = 0; n < d / 30; n++)
	{
		//  E = (6 n)!
		
		if (n == 0) E = Number(1).lvalue();
		
		else for (int j = 0; j < 6; j++)
		
		    E = E .multiply(6 * n - j);
		
		
		//  F = A + B n
		
		F = A .add(B.multiply(n));
		
		//  H[n] == n!^3 (3 n)! == (0!)^3 * 0!, (1!)^3 * 3!, (2!)^3 * 6!, ...
		//
		//       == H[n-1] n^3 ( 3 (n - 0) * 3(n - 1) * 3(n - 2) )
		
		if (n > 0)
		{
			H1 = H1 .divide(n * n) .divide(n);
			
			for (int j = 0; j < 3; j++)
			
			    H1 = H1 .divide(3 * n - j);
		}
		
		//  1/(12 pi) == (-1)^n E F G H ^-1
		
		product = (E) .multiply(F) .multiply(G) .multiply(H1);
		
		G = G .multiply(D2); //  G = C ^ (n+1/2)
		
		if ((n % 2) == 0)
		
		     sum = sum     .add(product);
		else sum = sum.subtract(product);
	}
	
	Number pi = sum .inverse() .divide(12);
	
	pi = pi .set_precision(digits);
	
	return *new Number(pi);
}


Number& Number::pow(double exp)
{
	Number number(exp);
	
	return this->pow(number);
}


Number& Number::pow(int exp)
{
	Number a (*this);
	
	Number y (1);
	
	int expsign = 0;
	
	if (exp < 0)
	
	    { exp *= -1; expsign = -1; }
	
	Number x = Number(exp) .abs();
	
	//  Compute y = a ^ x
	
	//  Use the square and multiply method
	
	while (!x.equals(0))
	{
		//  Accumulate squares
		
		if (x.test_bit(0))
		
		    y = a.multiply(y);
		
		//  Square the square
		
		a = a.square();
		
		//  Shift the exponent to the next bit
		
		x = x.shift_right(1);
	}
	
	if (expsign == -1) y = y.inverse();
	
	if ((this->sign == '-') && ((exp & 1) == 1))
	
	    y.sign = '-';
	
	return set_result(y);
}


Number& Number::pow(Number& exp)
{
	//  computes the function
	//
	//        x     (log a) x      k x
	//  y = a  ==  e          == e
	
	int precision = (this->precision >= exp.precision) ?
	                 this->precision :  exp.precision;
	
	if (exp.equals(0))
	
	    return set_result(Number(1).lvalue());
	
	if (exp.equals(1))
	
	    return set_result(*this);
	
	
	if (exp.is_complex())
	{
		//  e^(i z) == cos z + i sin z
		
		//  e^(a + ib) == e^a  e^(i b) == e^a (cos b + i sin b)
		
		Number a = exp.to_real();
		Number b = exp.to_imag();
		
		Number e2ib (exp.cos(), exp.sin());
		
		Number e2a = e_(precision).pow(a);
		
		Number product = e2a.multiply(e2ib);
		
		return set_result(product);
	}
	
	
	//  Compute the integer value
	
	Number y_int = this->pow1(exp.abs().to_integer());
	
	
	//  If exponent is an integer, return y_int
	
	if (exp.is_integer())
	
	    return set_result(y_int);
	
	
	//                                    k
	//  Solve for k in the equation a = e   or k = ln(a)
	
	Number base = Number(*this).set_precision(precision);
	
	Number k (Math::log(base.double_value()));
	
	
	//           k x                  k x        n
	//  Compute e   from the series  e   == (k x) / n!
	//
	//  where x is the exponent and n is the index 0, 1, 2, ...
	//
	//  This power series brings the fraction k x from the
	//  exponent down to the base to avoid doing square roots
	//  which are expensive.
	
	
	Number y_frac(0);
	
	Number factorial(1);
	
	
	//  iterations = digits / digits per iteration
	
	int d = 4 + (int) (2 * precision
	
	    / Math::log((double) precision, (double) 16));
	
	Number exp1 = exp.abs().to_fraction();
	
	Number kx = k.multiply(exp1) .set_precision(precision);
	
	Number kx2n = Number(1) .set_precision(precision);
	
	
	//  a^x == (k x) ^ n / n!
	
	for (int i = 1;  ; i++)
	{
		Number kx2ndf = kx2n .divide(factorial);
		
		y_frac = y_frac .add(kx2ndf);
		
		kx2n = kx2n .multiply(kx);
		
		factorial = factorial .multiply(i);
		
		if (kx2ndf.equals(0)) break;
	}
	
	//  Compute y = y_int * y_frac
	
	Number y = y_int.multiply(y_frac);
	
	
	//  Invert the output if the exponent is negative
	
	if (exp.less_than(0))
	
	    y = y.inverse();
	
	return set_result(y);
}


Number& Number::pow1(Number& exp)
{
	//  integer exponent method
	
	Number a (*this);
	
	Number y (1);
	
	Number x = Number(exp).abs();
	
	//  Use the square and multiply method
	
	while (!x.equals(0))
	{
		//  Accumulate squares
		
		if (x.test_bit(0))
		
		    y = a.multiply(y);
		
		//  Square the square
		
		a = a.square();
		
		//  Shift the exponent to the next bit
		
		x = x.shift_right(1);
	}
	
	if (exp.less_than(0))
	
	    y = y.inverse();
	
	if ((this->sign == '-') && exp.test_bit(0))
	
	    y.sign = '-';
	
	return set_result(y);
}



Number& Number::quad_divide(Number& divisor)
{

	//  The quadratic divider (or subtract and shift divider)
	//  uses the method of partial quotients to divide two numbers.
	
	//  This method is the inverse of the quadratic multiplier.
	//  Instead of doing small multiplies, adds and shifts, the
	//  quadratic divider does small divides, shifts and subtracts.
	//
	//  The quadratic divider works by dividing the most significant
	//  bits of the dividend by the most significant bits of the
	//  divisor, then multiplying the partial quotient by the divisor,
	//  shifting the partial quotient left or right until the most
	//  significant bits align with the dividend, subtracting the
	//  product from the dividend, and adding the partial quotient to
	//  the quotient until the dividend is reduced to a value less
	//  than the divisor.
	//
	//  Note that division (just like multiplication) can be done without
	//  a 32-bit multiplier or divider by shifting and subtracting bits.
	//  The processor's multiplier just makes the division operation
	//  several times faster because the dividend is reduced by several
	//  bits per iteration instead of only 1 to 2 bits per iteration.
	//  Doing multiplication or division by shifting and adding is
	//  equivalent to using a 1-bit multiplier.
	//
	//  To keep the quad_divide method simple, the method does not allow
	//  floating point dividends or divisors because the quad_divide method
	//  is called by the divide method which left shifts the dividend to
	//  the required precision (after equalizing the dividend and divisor
	//  lengths), removes the intpoints and precision, and then sets the
	//  quotient intpoint equal to the difference between the two intpoints.
	//  It would be redundant to replicate this code in the quad_divide method.
	
	
	Number dividend = *this;
	
	Number dividend1 = dividend;
	Number divisor1  = divisor;
	
	
	vector<int> vector0 = divisor1 .trim() .to_int_vector();
	
	long a, p = vector0[0];
	
	Number quotient(0);
	
	int counter = 1;
	
	
	while (true)
	{
		//  Divide the most significant int of the dividend
		//      by the most significant int of the divisor
		
		const vector<int> vector1 = dividend1 .trim() .to_int_vector();
		
		int vector_length1 = dividend1.vec_int.size();
		
		if (vector_length1 >= 2)
		
			a = (1L * vector1[0] & 0xffffffffL) * 0x100000000L
		          + (1L * vector1[1] & 0xffffffffL);
		else    a = (1L * vector1[0] & 0xffffffffL) * 0x100000000L;
		
		
		//  Maximize the number of significant digits in the quotient
		//  by left shifting the most significant digits in the dividend
		
		if      ((a & 0xffffffffffff0000L) == 0) a <<= 48;
		else if ((a & 0xffffffffff000000L) == 0) a <<= 40;
		else if ((a & 0xffffffff00000000L) == 0) a <<= 32;
		else if ((a & 0xffffff0000000000L) == 0) a <<= 24;
		else if ((a & 0xffff000000000000L) == 0) a <<= 16;
		else if ((a & 0xff00000000000000L) == 0) a <<=  8;
		
		if (a < 0)  a >>= 1;  a &= 0x7fffffffffffffffL;
		
		long q = a / p;
		
		while ((q & 0xffffffff00000000L) != 0)
		
		    { q >>= 2;  q &= 0x7fffffffffffffffL; }
		
		
		//  Compute the product of the partial quotient and the divisor
		
		Number product ( Math::multiply( divisor1.vec_int, (int) q ),
		
		    divisor.vec_int.size() + 1 );
		
		
		
		//  Subtract the product from the dividend if the product is > 0,
		//  or add the product to the dividend if the product is < 0, so that
		//  the most significant bits cancel out and the dividend contracts
		//  by several bits per iteration. (The sign of the dividend may
		//  alternate but the magnitude continues to decrease.) If the remainder
		//  or reduced dividend is negative when breaking out of the loop, then
		//  subtract 1 from the quotient because too much was subtracted from the
		//  dividend.
		
		//  Calculate the difference between the dividend and product lengths
		//  to left or right shift the product so that the most significant
		//  digits align and the dividend gets reduced by subtracting (or adding)
		//  the product.
		
		
		int d1 = (int) dividend1.bit_count() - 32*dividend.intpoint; // intpoint == 0
		int d2 = (int)   product.bit_count() - 32* product.intpoint; // intpoint == 0
		
		int d = d1 - d2;
		
		
		
		//  Q is a partial quotient
		
		Number Q (q);
		
		if      (d > 0) { Q = Q .shift_left(d, d); }
		else if (d < 0) { Q = Q .shift_right(-d) .trim(); }
		
		
		
		//  Q x divisor is a partial dividend
		
		//  Recalculate the product = Q x divisor
		//
		//  (with or without using multiplication)
		
		if (d > 0) product = product.shift_left(d, d);
		
		////  if (d < 0) product = Q .multiply(divisor1); // using multiplication
		if       (d < 0) product = product .shift_right((-d)
		              % (32*product.vec_int.size())) .trim();
		
		
		//  If the reduced dividend is less than the divisor, break out of the loop
		
		//  (Test d first to avoid the second test which is more expensive)
		
		if (d <= 0)
		{
			if (dividend1.equals(0)) break;
		}
		
		else if (d <  0)
		{
			if (dividend1.bit_count() <= divisor1.bit_count()) break;
		}
		
		
		//  Subtract the product from the dividend if the product is positive
		//  or add the product to the dividend if the product is negative
		//  so that the most significant bits cancel out and the dividend
		//  contracts by several bits per iteration. (The sign may alternate
		//  but the magnitude continues to decrease.)
		
		
		if (dividend1 .signum() == product .signum())
		{
			//  Subtract the (positive) product (once or twice) from the dividend
			
			dividend1 = dividend1 .subtract(product);
			
			quotient = quotient .add(Q);
			
			if (dividend1 .subtract(product) .bit_count() < dividend1.bit_count())
			{
				dividend1 = dividend1 .subtract(product);
				
				quotient = quotient .add(Q);
			}
		}
		
		else // if (dividend1 .signum() != product .signum())
		{
			//  Add the (negative) product (once or twice) to the dividend
			
			dividend1 = dividend1 .add(product);
			
			quotient = quotient .subtract(Q);
			
			if (dividend1 .add(product) .bit_count() < dividend1.bit_count())
			{
				dividend1 = dividend1 .add(product);
				
				quotient = quotient .subtract(Q);
			}
		}
		
		
		if (counter++ > dividend.bit_count() + 16)
		{
			string message = "quad divide method error";
			
			cout << message << endl; throw string(message);
		}
	}
	
	
	//  the dividend should contract around 16 bits per iteration except if the quotient is small
	
	long avg_bits_per_iter = (int) (dividend.bit_count() - dividend1.bit_count()) / counter;
	
	//  if (quotient.bit_count() > 64)  cout << "avg bits / iteration == " << avg_bits_per_iter;
	
	
	
	//  Correct for any small one-off errors
	
	while (this ->subtract(quotient.multiply(divisor)) .less_than(0))
	{
		quotient = quotient .subtract(1);  //  quotient is too large
	}
	
	Number remainder;
	
	do
	{	remainder = this ->subtract(quotient.multiply(divisor));
		
		//  Verify the addition, subtraction, and multiplication methods
		
		if (!remainder .add(divisor.multiply(quotient)) .equals(*this))
		{
			throw string("addition, subtraction, or multiplication error");
		}
		
		if (!remainder .less_than(divisor))
		{
			quotient = quotient .add(1);  //  quotient is too small
		}
	}
	
	while (!remainder .less_than(divisor));
	
	
	//  Verify the quotient
	
	remainder = this->subtract(quotient.multiply(divisor));
	
	if ( (remainder.less_than(0)) || (!remainder.less_than(divisor)) )
	{
		cout << "dividend - quotient x divisor == " << remainder.to_string(16) << endl;
		
		cout << "remainder compare zero == " << remainder.compare(0) << endl;
		
		cout << "remainder is less than zero == " <<
		
		    ((remainder.sign == '-') ? "true" : "false") << endl;
		
		cout << "dividend == " <<    this->to_string(16) << endl;
		cout << " divisor == " <<  divisor.to_string(16) << endl;
		cout << "quotient == " << quotient.to_string(16) << endl;
		
		string message = "quadratic division error";
		
		throw string(message);
	}
	
	return set_result(quotient);
}



Number& Number::root(int k)
{
	//  computes the kth root r of a number n using
	//
	//  Newton's iteration for root extraction
	
	//  Set r = number of int digits / k, then iterate
	//
	//  r == ((k-1) r + n / r ^ (k-1)) / k  or
	//
	//  r[i+1] == [(k-1) r[i] + n / r[i] ^ (k-1)] / k.
	//
	//  (The root r has to be set to the approximate size
	//  or else the iteration will converge at only one bit
	//  per iteration instead of converging quadratically
	//  or doubling the number of bits for each iteration.
	//
	//
	//  Example  Compute the square root r of a number n
	//
	//  (Newton's iteration for square root extraction)
	//
	//  Set k = 2, then iterate
	//
	//  r == [(k-1) r + n / r ^ (k-1)] / k
	//    == [(2-1) r + n / r ^ (2-1)] / 2
	//    ==      ( r + n / r ) / 2
	//
	//
	//  Example  Compute the minus one root r of a number n
	//
	//  (Newton's iteration for inversion)
	//
	//  Set k = -1, then iterate
	//
	//  r ==  [(k-1) r + n / r ^ (k-1)] / k
	//
	//    ==  [(-2) r + n / r ^ (-2)] / -1
	//
	//    ==  [ 2 r - n r ^ 2]
	//
	//    ==  r (2 - n r)
	//
	//
	//  Newton's iteration for k = -1 to 5
	//
	//  inverse root (n = r^-1)    r = (-2 r + n / r^-2) /-1 == r (2 - n r)
	//  zeroth  root (n = r^ 0)    r = (-1 r + n / r^-1) / 0
	//  linear  root (n = r^ 1)    r = ( 0 r + n / r^ 0) / 1
	//  square  root (n = r^ 2)    r = ( 1 r + n / r^ 1) / 2 == (r + n / r) / 2
	//  cube    root (n = r^ 3)    r = ( 2 r + n / r^ 2) / 3
	//  quartic root (n = r^ 4)    r = ( 3 r + n / r^ 3) / 4
	//  quintic root (n = r^ 5)    r = ( 4 r + n / r^ 4) / 5
	//  ...     ...                ...     ...
	//
	//  The inverse of v is computed by initializing u ~ 1 / v
	//  (so that u v ~ 1) and then iterating  u = u (2 - v u).
	//
	//  The square root of n is computed by initializing r ~ sqrt(n)
	//  (so that n ~ r^2) and then iterating  r = (r + n / r) / 2.
	
	
	if (this->equals(0) || this->equals(1))
	
	    //  zero and one are idempotent
	{
		return set_result(*this);
	}
	
	//  a ^  minus root == (1/a) ^ plus root
	//  a ^ zeroth root == a ^ (1/0) == infinity
	//  a ^  first root == a ^ (1/1) == a
	
	if (k <  0)  return this ->inverse() .root(-k);
	if (k == 0)  throw string("zeroth root");
	if (k == 1)
	{
		return set_result(*this);
	}
	
	
	//  Move the integer point so the number is close to 1
	//
	//  (or else the root will converge in O(n bits) instead of O(log(n bits))
	
	
	//  Remove the integer point from the number n = a / 2^b
	
	Number n =  this->trim();  n.intpoint = 0;
	
	int b = 32 * this->intpoint;
	
	
	//  Set the root r approx equal to n ^ (1/k)
	
	long bits = n.bit_count();
	
	Number r(1);  r = r .shift_left(bits/k, bits/k);
	
	
	//  The number of twos has to be divisible by k
	//
	//  so that the root of 2 ^ bits == 2 ^ (bits / k)
	
	while ((b % k) != 0)
	
	    { n = n .shift_left(1);  b++; }
	
	
	//  First iterate at 1/16 th of the precision
	
	int p = 8;
	
	for (int i = 0; i < 8; i++)
	{
		Number n_div_r2k_m1 = n.set_precision(p) .divide(r.pow(k-1));
		
		r = r .multiply(k-1) .add( n_div_r2k_m1 )
		
		    .divide(k) .trim();
	}
	
	
	//  Do a few more iterations doubling the precision for each iteration
	
	//  This is equivalent to 1.5 multiplications at full precision because
	//  all of the iterations except for the last iteration are done at much
	//  smaller precisions and these smaller multiplications are inexpensive
	//  because the sum of (1 / 2^k) ^ (1.58 or 2.00) == ... (1 / 8) ^ 1.58 +
	//  (1 / 4) ^ 1.58 + (1 / 2) ^ 1.58 + (1 / 1) ^ 1.58 < 1.5.
	
	//  Note that if the precision is set too high, it can make the method
	//  (and the computation of pi) several times slower
	
	for (   ; p < ((this->precision > 64) ? this->precision : 64); p *= 2)
	
	    r = r .multiply(k-1) .add( n.set_precision(p) .divide(r.pow(k-1)) ) .divide(k);
	
	
	//  the kth root of a (= n x 2^b) == root(n) x 2^(b/k)
	
	r = r .divide( Number(2) .pow(b/k) );
	
	
	//  Verify that the root is correct
	
	Number zero ("0", 16); zero = zero .set_precision(this->precision * 3/4);
	
	if (!r.pow(k) .subtract(*this) .abs() .equals(zero))
	{
		cout << "r^k == " << r.pow(k) .to_string(16);
		
		cout << "this == " << this->to_string(16);
		
		cout << "r^k / this == " << r.pow(k)
		
		    .divide(*this) .abs() .to_string(16) << " != 1";
		
		cout << "r^k - this == " << r.pow(k)
		
		    .subtract(*this) .abs() .to_string(16) << " != 0";
		
		cout << "r^k - this should have at least "
		
		    << zero.precision + " zero digits";
		
		string message = "root extraction error";
		
		throw string(message);
	}
	
	
	//  If the radicand is the power of an integer, then return the root r as an integer
	
	if (this->is_integer() && r.add(0.1).to_integer() .pow(k) .equals(*this))
	{
		Number root = r.add(0.1).to_integer();
		
		return set_result(root);
	}
	
	r = r.set_precision(n.precision);
	
	return set_result(r);
}


Number& Number::round()
{
	//  Round the number up or down
	
	Number number(*this);
	
	Number fraction = number.to_fraction();
	
	if (fraction.double_value() >= 0.5)
	{
		if (number.signum() >= 0)
		
		     number = number.to_integer()     .add(1);
		else number = number.to_integer().subtract(1);
	}
	
	number = number.to_integer();
	
	return set_result(number);
}


Number& Number::set_bit(long bit)
{
	//  The number is changed by the set_bit method
	
	//  and by the add_bit and clear_bit methods
	
	int size = this->vec_int.size();
	
	if ((bit / 32 + 1) > size)
	{
		//  Expand the vector
		
		int vec_length = (int) (bit) / 32 + 1;
		
		vector<int> vec_int (vec_length);
		
		for (int i = 0; i < vec_length; i++) vec_int[i] = 0;
		
		for (int i = 0; i < size; i++)
		
		    vec_int[vec_length - size] = this->vec_int[i];
		
		this->vec_int = vec_int;
	}
	
	Math::set_bit(this->vec_int, bit);
	
	return *this;
}




Number& Number::set_precision(int precision)
{

	//  sets the number of fractional digits and expands
	//  or truncates the right side of the number
	
	
	//  Example  2 int ints + 3 frac ints == 5 total ints
	//
	//  xxxxxxxxxxxxxxxx.xxxxxxxxxxxxxxxxxxxxxxxx
	//  |______________| |______________________|
	
	
	Number n (*this);
	
	n.precision = precision;
	
	int total_ints, frac_ints, int_ints;
	
	total_ints = n .vec_int.size();
	frac_ints = n .intpoint;
	int_ints = total_ints - frac_ints;
	
	if (int_ints < 0) int_ints = 0;
	
	if (int_ints == 0)  //  int value == 0
	{
		//  Expand the left side of the vector to equal 1 + intpoint
		
		int new_length = 1 + n.intpoint;
		
		vector<int> new_vec_int (new_length);
		
		for (int i = 0; i < new_length; i++) new_vec_int[i] = 0;
		
		//  Align the vectors and copy right (to left)
		
		for (int i = 0; i < n.vec_int.size(); i++)
		
		    new_vec_int[      new_length -1 -i] =
		      n.vec_int[n.vec_int.size() -1 -i];
		
		//  Re-assign the int vector
		
		n.vec_int = new_vec_int;
		
		total_ints = n.vec_int.size();
		 frac_ints = n.intpoint;
		  int_ints = total_ints - frac_ints;
		
		if (int_ints < 0) int_ints = 0;
	}
	
	
	if ((frac_ints < (precision + 7) / 8) && (precision > 0))
	{
		//  Expand the right side of the vector to equal the precision / 8
		
		int new_length = int_ints + 1 + precision / 8;
		
		vector<int> new_vec_int (new_length);
		
		for (int i = 0; i < new_length; i++) new_vec_int[i] = 0;
		
		//  Align the intpoints and copy left (to right)
		
		n.intpoint += new_length - n.vec_int.size();
		
		for (int i = 0; i < n.vec_int.size(); i++)
		
		    new_vec_int[i] = n.vec_int[i];
		
		//  Re-assign the vector
		
		n.vec_int = new_vec_int;
		
		total_ints = n.vec_int.size();
		 frac_ints = n.intpoint;
		  int_ints = total_ints - frac_ints;
		
		if (int_ints < 0) int_ints = 0;
	}
	
	
	else if (frac_ints > (precision + 7) / 8)
	{
		//  Truncate the right side of the vector to precision / 8
		
		int new_length = int_ints + 1 + precision / 8;
		
		vector<int> new_vec_int (new_length);
		
		for (int i = 0; i < new_length; i++) new_vec_int[i] = 0;
		
		int size = new_length <= n.vec_int.size() ?
			   new_length  : n.vec_int.size();
		
		//  Align the intpoints and copy left (to right)
		
		n.intpoint -= (n.vec_int.size() - new_length);
		
		for (int i = 0; i < size; i++)
		
		    new_vec_int[i] = n.vec_int[i];
		
		//  Re-assign the vector
		
		n.vec_int = new_vec_int;
		
		total_ints = n.vec_int.size();
		 frac_ints = n.intpoint;
		  int_ints = total_ints - frac_ints;
		
		if (int_ints < 0) int_ints = 0;
	}
	
	if (n .is_complex())
	{
		Number imag = n.to_imag();
		
		imag.set_precision(precision);
		
		n.vec_int1 = imag.vec_int1;
		n.precision1 = imag.precision1;
		n.sign1 = imag.sign1;
	}
	
	//  This test slows down the method significantly
	//  
	//  if (!this->to_integer() .equals(n.to_integer()))
	//  {
	//	string message = "set_precision error";
	//	
	//	cout << "this integer == " << this->to_integer().to_string(16) << endl;
	//	cout << "   n integer == " <<     n.to_integer().to_string(16) << endl;
	//	
	//	cout << message << endl; throw string(message);
	//  }
	
	return set_result(n);
}




Number& Number::shift_left(long bits)
{
	//  shifts left but does not expand
	//
	//  To expand and shift left, use shift_left(expansion, bits).
	
	Number n (*this);
	
	n.vec_int = Math::shift_left(n.vec_int, bits);
	
	return set_result(n);
}


Number& Number::shift_left(long expansion, long bits)
{
	//  expands and shifts left
	
	Number n (*this);
	
	n.vec_int = Math::shift_left(
	
	    n.vec_int, expansion, bits);
	
	return set_result(n);
}


Number& Number::shift_right(long bits)
{
	//  shifts right but does not contract the vector
	//
	//  To shift right and contract the vector, use
	//
	//  shift_right(bits) .trim() to remove the left zeros.
	
	Number n (*this);
	
	n.vec_int = Math::shift_right(n.vec_int, bits);
	
	return set_result(n);
}


int Number::signum()
{
	//  returns the sign of a number
	//
	//  as a number (1, 0, or -1)
	
	if (this->equals(0)) return 0;
	
	return (this->sign == '-') ? -1 : +1;
}


Number& Number::sin()
{
	return this->cos_sin()[1];
}

Number& Number::sinh()
{
	//  returns the hyperbolic sine
	
	//  The hyperbolic sine is defined by
	//                         u     -u
	//  sinh(u)  =  1 / 2  ( e   -  e  )
	//
	//  The hyperbolic sine is related to the sine by
	//
	//  i sin(u)  == sinh(i u)  and  i sinh(u) == sin(i u)
	
	
	this->precision = (this->precision > 0) ? this->precision : 16;
	
	Number sinh = e_(this->precision).pow(*this)
	         .add(e_(this->precision).pow( this->negate())).divide(2);
	
	return set_result(sinh);
}

Number& Number::sqrt()
{
	//  returns the square root
	
	return this->root(2);
}

Number& Number::square()
{
	//  squares the number
	
	Number n (*this);
	
	Number square = n .multiply(n);
	
	return set_result(square);
}

Number& Number::subtract(int subtrahend)
{
	Number number(subtrahend);
	
	return this->subtract(number);
}

Number& Number::subtract(long subtrahend)
{
	Number number(subtrahend);
	
	return this->subtract(number);
}

Number& Number::subtract(double subtrahend)
{
	Number number(subtrahend);
	
	return this->subtract(number);
}


Number& Number::subtract(Number& subtrahend)
{
	//  subtracts two numbers
	
	Number number(subtrahend);
	
	return this->add(number.negate());
}


Number& Number::tan()
{

	//  computes the tangent of x
	
	//  tan(x) = sin(x) / cos(x)
	
	int p = this->precision;
	
	Number zero = Number(0).set_precision(p);
	
	if (this->cos() .equals(zero))
	{
		string message = "tan(pi/2) is undefined";
		
		cout << message << endl; throw string(message);
	}
	
	//  Compute the tangent of this
	
	return this->sin() .divide(this->cos());
}


bool Number::test_bit(long bit)
{
	if ((1 + bit / 32) > this->length())
	
	    return false;
	
	return Math::test_bit(this->vec_int, bit);
}


string Number::to_alphabetical_string()
{

	//  converts a number to an alphabetical string
	
	
	//  Example  Convert the number 314159.26 to an alphabetical string
	//
	//  cout << Number("314159.26") .to_alphabetical_string() << endl;
	//
	//  three hundred fourteen thousand, one hundred fifty-nine
	//
	//
	//  Example  Convert the number
	//
	//  123456789012345678901234567890123456789012345678901234567890123456 to an alphabetical string
	//
	//  cout << Number("123456789012345678901234567890123456789012345678901234567890123456"
	//
	//  	    .to_alphabetical_string() << endl;
	//
	//  One hundred twenty-three vigintillion, four hundred fifty-six novemdecillion,
	//  seven hundred eighty-nine octodecillion, twelve septendecillion, three hundred
	//  forty-five sexdecillion, six hundred seventy-eight quindecillion, nine hundred
	//  one quattuordecillion, two hundred thirty-four tredecillion, five hundred
	//  sixty-seven duodecillion, eight hundred ninety undecillion, one hundred
	//  twenty-three decillion, four hundred fifty-six nonillion, seven hundred
	//  eighty-nine octillion, twelve septillion, three hundred forty-five sextillion,
	//  six hundred seventy-eight quintillion, nine hundred one quadrillion, two hundred
	//  thirty-four trillion, five hundred sixty-seven billion, eight hundred ninety million,
	//  one hundred twenty-three thousand, four hundred fifty-six
	
	
	string base0[10] =
	
	    {" ", "one ", "two ", "three ", "four ", "five ", "six ", "seven ", "eight ", "nine "};
	
	string base1[10] =
	
	    {"ten ", "eleven ", "twelve ", "thirteen ", "fourteen ", "fifteen ", "sixteen ",
	
		"seventeen ", "eighteen ", "nineteen "};
	
	string base2[10] = {" ", " ", "twenty", "thirty", "forty", "fifty", "sixty", "seventy",
	
	    "eighty", "ninety"};
	
	string _hundred = "hundred ";
	
	
	string base3[22] = {"", "thousand", "million", "billion", "trillion", "quadrillion",
	
	    "quintillion", "sextillion", "septillion", "octillion", "nonillion", "decillion",
	    "undecillion", "duodecillion", "tredecillion", "quattuordecillion", "quindecillion",
	    "sexdecillion", "septendecillion", "octodecillion", "novemdecillion", "vigintillion"};
	
	
	
	//  If the number is greater than 999.999... vigintillion
	//
	//  	return the empty string
	
	string max_val_str =
	
	    string("999999999999999999999999999999999")
	  + string("999999999999999999999999999999999");
	
	Number max_val (max_val_str);
	
	if (this->greater_than(max_val)) return "";
	
	string n = this ->add(0.1)
	
	    .to_integer() .to_string(10);
	
	while ((n.length() % 3) != 0)
	
	    n = string("0") + n;
	
	int array_length = n.length() / 3;
	
	string   group[array_length];
	string sbarray[array_length];
	
	for (int i = 0; i < array_length; i++)
	{
		group[i] = n.substr(3 * i, 3 * (i + 1) - (3 * i));
		
		sbarray[i] = string("");
		
		char c = group[i][0];
		
		if (c != '0')
		{
			sbarray[i].append(base0[c - '0']);
			sbarray[i].append(_hundred);
		}
		
		c = group[i][1];
		
		if (c > '1')
		{
			sbarray[i].append(base2[c - '0']);
			
			if ((group[i][2] - '0') == 0)
			
			     sbarray[i].append(string(" ") + base0[group[i][2] - '0']);
			else sbarray[i].append(string("-") + base0[group[i][2] - '0']);
		}
		
		else if (c == '1') sbarray[i].append(base1[group[i][2] - '0']);
		else if (c == '0') sbarray[i].append(base0[group[i][2] - '0']);
		
		group[i] = string(sbarray[i]);
	}
	
	string sb = string("");
	
	for (int i = 0; i < array_length; i++)
	{
		if (group[i].find_first_not_of(" ") != string::npos)
		{
			if (i != 0) sb.append(", ");
			
			group[i] += base3[(array_length -i -1)];
			
			sb.append(group[i]);
		}
	}
	
	if ((sb.length() > 0) && (sb[0] >= 'a') && (sb[0] <= 'z'))
	
	    sb[0] = (char) (sb[0] - 'a' + 'A');
	
	string str = string(sb);
	
	
	int pos = 0;
	
	//  Replace double spaces with single spaces
	
	while ((pos = str.find("  ")) != string::npos)
	
	    str.erase(pos, 1);
	
	return str;
}


Number& Number::to_fraction()
{
	//  returns the fractional part
	
	return this->subtract(this->to_integer());
}


Number& Number::to_imag()
{
	//  returns the imag component
	//  of a complex number
	
	if (!this->is_complex())
	
	    return set_result(Number(0).lvalue());
	
	Number number (this->vec_int1);
	
	number.intpoint  = this->intpoint1;
	number.precision = this->precision1;
	number.sign      = this->sign1;
	
	return set_result(number);
}


vector<int> Number::to_int_vector()
{
	//  returns a copy of the int vector of a number
	
	//  (The int vector is just the magnitude of the number.
	//  It does not include the intpoint, precision, or sign.
	//  Methods that use twos complement arithmetic can
	//  interpret the sign from the zeroth int.)
	
	return this->vec_int;
}


vector<int> Number::to_int_vector(int size)
{
	//  returns a padded copy of the int vector of a number
	
	//  (If size > vector size, then the vector is padded w/zeros)
	
	int ints = (size >= this->vec_int.size()) ?
	
	    size : this->vec_int.size();
	
	vector<int> vec_int (ints);
	
	for (int i = 0; i < vec_int.size(); i++) vec_int[i] = 0;
	
	int offset = ints - this->vec_int.size();
	
	for (int i = 0; i < this->vec_int.size(); i++)
	
	    vec_int[offset + i] = this->vec_int[i];
	
	return vec_int;
}


Number& Number::to_integer()
{
	
	if ((this->precision == 0) && (this->intpoint == 0))
	{
		return this->set_result(*this);
	}
	
	int vec_length = this->length() - this->intpoint;
	
	if (vec_length <= 0)
	{
		return set_result(Number(0).lvalue());
	}
	
	Number n (Math::shift_right(vec_int, 32 * intpoint));
	
	n = n .trim();
	
	n.sign = this->sign;
	
	return set_result(n);
}


Number& Number::to_real()
{
	//  returns the real component
	//  of a complex number
	
	Number n (this->vec_int);
	
	n.intpoint  = this->intpoint;
	n.precision = this->precision;
	n.sign      = this->sign;
	
	return set_result(n);
}


string Number::to_string()
{
	return this->to_string(0, 10);
}

string Number::to_string(int radix)
{
	return this->to_string(0, radix);
}

string Number::to_string(int digits, int radix)
{
	//  converts a number to string in any base or radix <= 16
	
	if (this->is_complex())
	{
		//  Return the complex number string
		
		Number real = this->to_real();
		Number imag = this->to_imag();
		
		string rstring, istring;
		
		rstring = real.to_string(digits, radix);
		istring = imag.to_string(digits, radix);
		
		if (imag.sign == '+')
		
		     return rstring + string(" + ") + istring + string(" i");
		else return rstring + string(" - ") + istring + string(" i");
	}
	
	
	int maxsize = 256; //  256 ints * 8 digits / int == 2048 digits
	
	//  If the number of digits is greater than maxsize, this method
	//
	//  will call another to_string method to convert the number to string
	
	
	//  Make a copy the number
	
	Number number(*this);
	
	
	//  Expand the vector if necessary to make the leading zeros explicit
	
	if (number.intpoint > number.vec_int.size())
	
	    number = number .set_precision(this->precision);
	
	
	//  If base is > 16, throw an illegal argument exception
	
	if ((radix < 2) || (radix > 16))
	{
		//  The to_string method may only be accurate to base 27;
		//
		//  it was only designed for base 2 to 16
		
		string message = "Base must be between 2 and 16";
		
		message += string("\n\nuse the to_int_vector(long radix) method for bases > 36");
		
		cout << message << endl; throw string(message);
	}
	
	string  integerstring = "0";
	string fractionstring = "0";
	
	string str = "0";  int index = 0;
	
	
	vector<int> vec_int = number.vec_int; // 32-bit ints
	
	int vec_length = number.vec_int.size();
	
	vector<char> char_256_vector = Convert::int_vector_to_char_vector(vec_int);
	
	int char_256_length = 4 * vec_length; // 8-bit chars
	
	vector<char> char_16_vector = Convert::char_256_vector_to_char_16_vector(char_256_vector);
	
	int char_16_length = 2 * char_256_length; // 4-bit chars
	
	for (int i = 0; i < char_16_length; i++)
	
	    str += char_16_vector[i];
	
	
	
	if (radix == 16)
	{
	
		integerstring = str.substr(0, str.length() - 8 * number.intpoint);
		
		if ((integerstring == "") || integerstring.empty()) integerstring = "0";
		
		fractionstring = str.substr(str.length() - 8 * number.intpoint, str.length());
		
		
		//  Set the precision
		
		if ((this->precision > 0) && (this->precision < fractionstring.length()))
		
		    fractionstring = fractionstring.substr(0, this->precision);
		
		
		//  Remove fraction string zeros
		
		index = 0;
		
		while ((fractionstring.length() -1 - index >= 0)
		
		    && (fractionstring[fractionstring.length() -1 - index] == '0')
		
			&& (fractionstring.length() - index > number.precision))
		{
			index += 1;
		}
		
		fractionstring = fractionstring.substr(0, fractionstring.length() - index);
		
		
		//  Add fraction string zeros
		
		index = this->precision - fractionstring.length();
		
		string sb = string(fractionstring);
		
		while ((index / 8) > 0)
		{
			sb.append("00000000");
			
			index -= 8;
		}
		
		while ((index / 1) > 0)
		{
			sb.append("0");
			
			index -= 1;
		}
		
		fractionstring = string(sb);
	}
	
	
	
	else if (radix != 16)
	{
	
		if ((str.length() - number.intpoint * 8) <= 0)
		
		     integerstring = "0";
		
		else integerstring = str.substr(0, str.length() - number.intpoint * 8);
		
		
		if (this->intpoint != 0)
		{
			fractionstring = str.substr(
			
			    str.length() - number.intpoint * 8, str.length());
			
			string fractionstring1(fractionstring);
		}
		
		
		//  Convert radix-16 integer or fraction string to radix-n string
		
		if (!integerstring.empty() && !Number(integerstring, 16).equals(0))
		{
			while (true)
			{
				//  2^k & (2^k -1) == 10000... & 01111... == 0
				
				if ((radix & (radix - 1)) != 0)  // radix != power of 2
				
				    integerstring = Number(integerstring, 16) .add(1) .to_string(16);
				
				if (integerstring.length() / 8 > maxsize)
				{
					//  Use recursion to convert a large hex string to base radix
					
					integerstring = Number(integerstring, 16)
					
					    .subtract(1) .to_string1(radix);
					
					break;
				}
				
				int length = 4 + (int) (integerstring.length()
				
				    * (Math::log(16.0) / Math::log(radix)));
				
				number = Number(integerstring, 16) .lvalue();
				
				number = number .set_precision(length);
				
				Number base (radix);
				
				Number base_to_length = base .pow(length);
				
				number = number .divide( base_to_length );
				
				string integerstring1 = "";
				
				
				//  Remove zeros
				
				while (number.multiply(radix).int_value() == 0)
				{
					number = number.multiply(radix);
					
					length = length - 1;
				}
				
				for (int i = 0; i < length; i++)
				{
					number = number.multiply(radix);
					
					int digit = number.int_value();
					
					if (digit < 10)
					
					     integerstring1 += (char) (digit -  0 + '0');
					else integerstring1 += (char) (digit - 10 + 'a');
					
					number = number.subtract(digit);
				}
				
				
				//  Verify that the integer string is correct
				
				Number number1 (integerstring1, radix);
				Number number2 (integerstring, 16);
				
				//  If radix is not a power of 2, add 1 before comparing
				
				//  2^k & (2^k -1) == 10000... & 01111... == 0
				
				if ((radix & (radix - 1)) != 0)
				
				    number1 = number1 .add(1);
				
				if (!number1.equals(number2))
				{
					string message = "error in to_string() method";
					
					message += string("\nnumber1 = ") + number1 .to_string(16);
					message += string("\nnumber2 = ") + number2 .to_string(16);
					
					cout << message << endl; throw string(message);
				}
				
				integerstring = integerstring1;
				
				break;
			}
		}
		
		
		if (!fractionstring.empty() && !Number(fractionstring, 16).equals(0))
		{
		
			while (true)
			{
				int length = fractionstring.length();
				
				if (fractionstring.length() / 8 > maxsize)
				{
					//  Use recursion to convert the string
					
					fractionstring = Number(fractionstring, 16)
					
					    .to_string2(length, radix);
					
					break;
				}
				
				number = Number(fractionstring, 16).set_precision(length)
				
				    .divide(Number(16).pow(length));
				
				fractionstring = "";
				
				for (int i = 0; i < this->precision + 4; i++)
				{
					number = number.multiply(radix);
					
					int digit = number.int_value();
					
					if (digit < 10)
					
					     fractionstring += (char) (digit -  0 + '0');
					else fractionstring += (char) (digit - 10 + 'a');
					
					number = number.subtract(digit);
				}
				
				break;
			}
		}
		
		number.sign = this->sign;
	}
	
	
	//  Set the length of the fraction string equal to the precision
	
	if ((this->precision > 0) && (this->precision + 1 < fractionstring.length()))
	
	    fractionstring = fractionstring.substr(0, this->precision + 1);
	
	
	//  Remove leading zeros
	
	//  ...
	
	//  ...
	
	
	int pos = 0;
	
	while ((pos = str.find(" ")) != string::npos)
	
	    str .erase(pos, 1);
	
	if (str.empty()) str = "0";
	
	
	
	//  Join the integer and fraction strings
	
	str = (integerstring.empty() ? "0" : integerstring)
	
	    + (((this->precision == 0) || fractionstring.empty()) ?
	
		"" : string(".") + fractionstring);
	
	//  If string starts with a period, prepend a zero
	
	if (str.find_first_of(".") == 0) str = string("0") + str;
	
	
	//  Pad the left side of the integer to the minimum number of digits
	
	//  If the radix is 16, the number will be padded with zeros because
	//  base 16 is read by computers and is used for cryptography. If the
	//  radix is not a power of 2, the number will be padded with spaces.
	
	int length = str.length();
	
	if (digits > length)
	{
		index = 0;  string sb;
		
		while (digits - index++ > length)
		{
			if ((radix == 16) || (radix == 8) || (radix == 4) || (radix == 2))
			
			     sb.append("0");
			else sb.append(" ");
		}
		
		sb.append(str);
		
		str = string(sb);
	}
	
	
	//  Pad the right side of the number to precision + 1 digits
	
	length = fractionstring.length();  index = 0;
	
	if ((this->precision > 0) && (this->precision + 1 > length))
	
	    while (this->precision + 1 - index++ > length)
	
		str.append("0");
	
	
	//  Round the fraction up one digit
	
	if ((radix <= 16) && (this->precision > 0) && (this->vec_int[this->length() -1] != 0))
	{
		string sb = string(str);
		
		sb.insert(0, " "); // space for carry bit
		
		length = sb.length();  bool carry = false;
		
		for (int i = length - 1; i >= 1; i--)
		{
			int c = sb[i];
			
			if (c == '.')  continue;
			
			//  Roll the digit over to zero or add one and break
			
			if ( ((c - '0' +  0) == (radix -1))
			  || ((c - 'a' + 10) == (radix -1)) )
			
			     { sb[i] = '0';  carry = true; }
			
			else { sb[i] = (char) (c + 1);  carry = false; }
			
			if (carry == false) break;
		}
		
		if (carry == true)  sb[0] = '1';
		
		int pos = 0;
		
		while (sb.find(" ") != string::npos)
		
		    sb .erase(pos, 1);
		
		str = string(sb);
	}
	
	
	//  Remove the last fraction digit
	
	if (this->precision > 0)
	
	    str = str.substr(0, str.length() - 1);
	
	
	//  Prepend the sign
	
	//  By convention, the plus sign is omitted
	
	if (this->sign == '+') ;
	
	////    str = " " + str;
	
	else if (this->sign == '-')
	{
		//  Zero can be plus or minus, but no sign is displayed for zero
		
		str = string("-") + str;
	}
	
	
	//  Verify that | Number(string, radix) - this | < 1
	
	//  ...   ...
	
	//  ...   ...
	
	
	//  Return the number string
	
	return str;
}



Number& Number::trim()
{

	//  trims the leading zero ints
	//
	//  If the number array value equals zero,
	//  the vector will be trimmed to one int.
	
	
	//  Trim the leading zeros
	
	Number n (*this);
	
	n.vec_int = Math::trim(n.vec_int);
	
	
	if (!this->is_complex())
	
	    return set_result(n);
	
	
	//  else  complex
	{
		//  ...
	}
}



























double Math::abs(double d)
{
	return d >= 0 ? d : d*-1.0;
}

int Math::abs(int i)
{
	return i >= 0 ? i : i*-1;
}

long Math::abs(long l)
{
	return l >= 0 ? l : l*-1L;
}


vector<int> Math::add(vector<int> addend, vector<int> summand)
{
	//  adds two signed vectors
	
	int  addend_length =  addend.size();
	int summand_length = summand.size();
	
	int length = addend_length > summand_length ? addend_length : summand_length;
	
	int a_array[length],  b_array[length];
	
	//  Initialize the vectors with the sign int so that the
	//  sign of the smaller vector will expand to the left
	
	int a_sign_int = ( addend[0] < 0) ? -1 : 0;
	int b_sign_int = (summand[0] < 0) ? -1 : 0;
	
	for (int i = 0; i < length; i++) a_array[i] = a_sign_int;
	for (int i = 0; i < length; i++) b_array[i] = b_sign_int;
	
	//  Copy right (to left)
	
	for (int i = 0; i < addend_length; i++)
	
	    a_array[length -1 -i] = addend[addend_length -1 -i];
	
	for (int i = 0; i < summand_length; i++)
	
	    b_array[length -1 -i] = summand[summand_length -1 -i];
	
	int *c_array = add(a_array, b_array, length);
	
	//  Number a_test (Convert::int_array_to_int_vector(a_array, length), true);
	//  Number b_test (Convert::int_array_to_int_vector(b_array, length), true);
	//  Number c_test (Convert::int_array_to_int_vector(c_array, length), true);
	
	vector<int> vec_int (length);
	
	for (int i = 0; i < length; i++)
	
	    vec_int[i] = c_array[i];
	
	delete[] c_array;
	
	return vec_int;
}



int *Math::add(int addend[], int summand[], int elements)
{

	//  adds two int arrays and returns an int array
	//
	//  This method is used by the Number class
	
	
	int a_length = elements;
	int b_length = elements;
	int c_length = elements;
	
	int a[elements], b[elements];
	
	int *c = new int[elements];
	
	for (int i = 0; i < a_length; i++) a[i] =  addend[i];
	for (int i = 0; i < b_length; i++) b[i] = summand[i];
	for (int i = 0; i < c_length; i++) c[i] = 0;
	
	
	long temp, temp_a, temp_b;
	
	temp = temp_a = temp_b = 0L;
	
	for (int i = 0; i < elements; i++)
	{
		temp_a += a[a_length -1 -i];
		temp_b += b[b_length -1 -i];
		
		temp_a &= 0xffffffffL;
		temp_b &= 0xffffffffL;
		
		temp += temp_a;
		temp += temp_b;
		
		c[c_length -1 -i] += temp;
		
		temp >>= 32;
		
		temp_a = temp_b = 0L;
	}
	
	return c;
}


void Math::add_bit(vector<int>& vec_int, long bit)
{
	//  adds a bit to a vector
	
	long n = bit;
	
	int elements = vec_int.size();
	
	if ((bit / 32) >= 1)
	{
		//  Roll over the -1 ints
		
		while (((elements -1 - n / 32) == -1)
		
		   && (vec_int[elements -1 - n / 32] == -1))
		     { vec_int[elements -1 - n / 32] =   0; n += 32; }
	}
	
	while ((n < 32 * elements)
	
	   //  Roll over the 1 (or "-1") bits
	
	   && (get_bit(vec_int, n) == 1))
	     clear_bit(vec_int, n++);
	
	
	if (n < 32 * elements)
	
	    set_bit(vec_int, n);
}



void Math::add_bit(int *array, int elements, long bit)
{
	//  adds a bit to an int array
	
	long n = bit;
	
	if ((bit / 32) >= 1)
	{
		//  Roll over the -1 ints
		
		while (((elements -1 - n / 32) == -1)
		
		   && (array[elements -1 - n / 32] == -1))
		     { array[elements -1 - n / 32] =   0; n += 32; }
	}
	
	while ((n < 32 * elements)
	
	   //  Roll over the 1 (or "-1") bits
	
	   && (get_bit(array, elements, n) == 1))
	     clear_bit(array, elements, n++);
	
	if (n < 32 * elements)
	
	    set_bit(array, elements, n);
}


vector<int> Math::and_(vector <int> a, vector<int> b)
{
	//  multiplies bits (modulo 2)
	
	int a_length = a.size();
	int b_length = b.size();
	
	int size = a_length < b_length ?
                   a_length : b_length;
	
	vector<int> c (size);
	
	for (int i = 0; i < size; i++) c[i] = 0;
	
	for (int i = 0; i < size; i++)
	
	    c [size -1 -i] = a[a_length -1 -i] & b[b_length -1 -i];
	
	return  c;
}


int *Math::and_(int a[], int b[], int elements1, int elements2)
{
	//  multiplies bits (modulo 2)
	
	int a_length = elements1;
	int b_length = elements2;
	
	int size = a_length < b_length ?
                   a_length : b_length;
	
	int * c = new int[size];
	
	for (int i = 0; i < size; i++) c[i] = 0;
	
	for (int i = 0; i < size; i++)
	
	    c [size -1 -i] = a[a_length -1 -i] & b[b_length -1 -i];
	
	return  c;
}


int Math::bit_count(int n)
{
	if (n < 0) return 32;
	
	int j = 0;
	
	while (n != 0)
	
	    { n >>= 1; j++; }
	
	return j;
}


int Math::bit_count(long n)
{
	if (n < 0) return 64;
	
	int j = 0;
	
	while (n != 0)
	
	    { n >>= 1; j++; }
	
	return j;
}



long Math::bit_count(vector<int> vec_int)
{
	int i = 0, j = 0;
	
	int vec_length = vec_int.size();
	
	//  Count the leading zero ints
	
	while ((i < vec_length - 1) && (vec_int[i] == 0)) i++;
	
	int b = vec_int[i];
	
	if (b < 0) j = 32;
	
	else while (b != 0) { b >>= 1; j++; }
	
	int bits = 32*(vec_length -1 -i) + j;
	
	return bits;
}



long Math::bit_count(int array[], int elements)
{
	int i = 0, j = 0;
	
	int array_length = elements;
	
	//  Count the leading zero ints
	
	while ((i < array_length - 1) && (array[i] == 0)) i++;
	
	int b = array[i];
	
	if (b < 0) j = 32;
	
	else while (b != 0) { b >>= 1; j++; }
	
	int bits = 32*(array_length -1 -i) + j;
	
	return bits;
}


void Math::clear_bit(vector<int>& vec_int, long bit)
{
	unsigned int bit_mask[32] =
	{
	      0xfffffffe, 0xfffffffd, 0xfffffffb, 0xfffffff7,
	      0xffffffef, 0xffffffdf, 0xffffffbf, 0xffffff7f,
	      0xfffffeff, 0xfffffdff, 0xfffffbff, 0xfffff7ff,
	      0xffffefff, 0xffffdfff, 0xffffbfff, 0xffff7fff,
	      0xfffeffff, 0xfffdffff, 0xfffbffff, 0xfff7ffff,
	      0xffefffff, 0xffdfffff, 0xffbfffff, 0xff7fffff,
	      0xfeffffff, 0xfdffffff, 0xfbffffff, 0xf7ffffff,
	      0xefffffff, 0xdfffffff, 0xbfffffff, 0x7fffffff,
	};
	
	vec_int[vec_int.size() -1 - (int) (bit / 32)] &= bit_mask[bit % 32];
}



void Math::clear_bit(int array[], int elements, long bit)
{
	unsigned int bit_mask[32] =
	{
	      0xfffffffe, 0xfffffffd, 0xfffffffb, 0xfffffff7,
	      0xffffffef, 0xffffffdf, 0xffffffbf, 0xffffff7f,
	      0xfffffeff, 0xfffffdff, 0xfffffbff, 0xfffff7ff,
	      0xffffefff, 0xffffdfff, 0xffffbfff, 0xffff7fff,
	      0xfffeffff, 0xfffdffff, 0xfffbffff, 0xfff7ffff,
	      0xffefffff, 0xffdfffff, 0xffbfffff, 0xff7fffff,
	      0xfeffffff, 0xfdffffff, 0xfbffffff, 0xf7ffffff,
	      0xefffffff, 0xdfffffff, 0xbfffffff, 0x7fffffff,
	};
	
	array[elements -1 - (int) (bit / 32)] &= bit_mask[bit % 32];
}



int Math::compare(vector<int> a, vector<int> b)
{

	//  Unsigned int comparator  returns 0, 1, or -1
	
	//  Note that the java and C++ comparators are slightly different
	//  because in C++  an int & 0x80000000 == (a signed int) & (an unsigned int),
	//  whereas in java an int & 0x80000000 == (a signed int) & (   a signed int).
	//
	//  In java the result of the && equals a negative value whereas in java it
	//  equals a positive value. The two versions would be identical if the C++
	//  version used (signed) 0x80000000.
	
	int a_length = a.size();
	int b_length = b.size();
	
	int d = abs((int) (a_length - b_length));
	
	if (a_length > b_length)
	{
		for (int i = 0; i < d; i++) if (a[i] != 0) return 1;
		
		if      ((a[d] & 0x80000000) > (b[0] & 0x80000000)) return  1;
		else if ((a[d] & 0x80000000) < (b[0] & 0x80000000)) return -1;
		else if ((a[d] & 0x7fffffff) > (b[0] & 0x7fffffff)) return  1;
		else if ((a[d] & 0x7fffffff) < (b[0] & 0x7fffffff)) return -1;
		
		else for (int i = 0; i < b_length; i++)
		{
			if      ((a[d+i] & 0x80000000) > (b[i] & 0x80000000)) return  1;
			else if ((a[d+i] & 0x80000000) < (b[i] & 0x80000000)) return -1;
			else if ((a[d+i] & 0x7fffffff) > (b[i] & 0x7fffffff)) return  1;
			else if ((a[d+i] & 0x7fffffff) < (b[i] & 0x7fffffff)) return -1;
	    	}
	}
	
	
	else if (b_length > a_length)
	{
		for (int i = 0; i < d; i++) if (b[i] != 0) return -1;
		
		if      ((a[0] & 0x80000000) > (b[d] & 0x80000000)) return  1;
		else if ((a[0] & 0x80000000) < (b[d] & 0x80000000)) return -1;
		else if ((a[0] & 0x7fffffff) > (b[d] & 0x7fffffff)) return  1;
		else if ((a[0] & 0x7fffffff) < (b[d] & 0x7fffffff)) return -1;
		
		else for (int i = 0; i < a_length; i++)
		{
			if      ((a[i] & 0x80000000) > (b[d+i] & 0x80000000)) return  1;
			else if ((a[i] & 0x80000000) < (b[d+i] & 0x80000000)) return -1;
			else if ((a[i] & 0x7fffffff) > (b[d+i] & 0x7fffffff)) return  1;
			else if ((a[i] & 0x7fffffff) < (b[d+i] & 0x7fffffff)) return -1;
		}
	}
	
	else if (a_length == b_length)
	{
		if      ((a[0] & 0x80000000) > (b[0] & 0x80000000)) return  1;
		else if ((a[0] & 0x80000000) < (b[0] & 0x80000000)) return -1;
		else if ((a[0] & 0x7fffffff) > (b[0] & 0x7fffffff)) return  1;
		else if ((a[0] & 0x7fffffff) < (b[0] & 0x7fffffff)) return -1;
		
		else for (int i = 0; i < a_length; i++)
		{
			if      ((a[i] & 0x80000000) > (b[i] & 0x80000000)) return  1;
			else if ((a[i] & 0x80000000) < (b[i] & 0x80000000)) return -1;
			else if ((a[i] & 0x7fffffff) > (b[i] & 0x7fffffff)) return  1;
			else if ((a[i] & 0x7fffffff) < (b[i] & 0x7fffffff)) return -1;
		}
	}
	
	return 0;
}



int Math::compare(int a[], int b[], int elements1, int elements2)
{
	//  Unsigned int comparator  returns 0, 1, or -1
	
	int a_length = elements1;
	int b_length = elements2;
	
	int d = abs((int) (a_length - b_length));
	
	if (a_length > b_length)
	{
		for (int i = 0; i < d; i++) if (a[i] != 0) return 1;
		
		if      (((a[d] & 0x80000000) - (b[0] & 0x80000000)) < 0)  return  1;
		else if (((a[d] & 0x80000000) - (b[0] & 0x80000000)) > 0)  return -1;
		else if (((a[d] & 0x7fffffff) - (b[0] & 0x7fffffff)) > 0)  return  1;
		else if (((a[d] & 0x7fffffff) - (b[0] & 0x7fffffff)) < 0)  return -1;
		
		else
		{	for (int i = 0; i < b_length; i++)
			{
				if      (((a[d + i] & 0x80000000) - (b[i] & 0x80000000)) < 0)  return  1;
				else if (((a[d + i] & 0x80000000) - (b[i] & 0x80000000)) > 0)  return -1;
				else if (((a[d + i] & 0x7fffffff) - (b[i] & 0x7fffffff)) > 0)  return  1;
				else if (((a[d + i] & 0x7fffffff) - (b[i] & 0x7fffffff)) < 0)  return -1;
			}
		}
	}
	
	else if (b_length > a_length)
	{
		for (int i = 0; i < d; i++) if (b[i] != 0) return -1;
		
		if      (((a[0] & 0x80000000) - (b[d] & 0x80000000)) < 0)  return  1;
		else if (((a[0] & 0x80000000) - (b[d] & 0x80000000)) > 0)  return -1;
		else if (((a[0] & 0x7fffffff) - (b[d] & 0x7fffffff)) > 0)  return  1;
		else if (((a[0] & 0x7fffffff) - (b[d] & 0x7fffffff)) < 0)  return -1;
		
		else
		{	for (int i = 0; i < a_length; i++)
			{
				if      (((a[i] & 0x80000000) - (b[d + i] & 0x80000000)) < 0)  return  1;
				else if (((a[i] & 0x80000000) - (b[d + i] & 0x80000000)) > 0)  return -1;
				else if (((a[i] & 0x7fffffff) - (b[d + i] & 0x7fffffff)) > 0)  return  1;
				else if (((a[i] & 0x7fffffff) - (b[d + i] & 0x7fffffff)) < 0)  return -1;
			}
		}
	}
	
	else if (a_length == b_length)
	{
		if      (((a[0] & 0x80000000) - (b[0] & 0x80000000)) < 0)  return  1;
		else if (((a[0] & 0x80000000) - (b[0] & 0x80000000)) > 0)  return -1;
		else if (((a[0] & 0x7fffffff) - (b[0] & 0x7fffffff)) > 0)  return  1;
		else if (((a[0] & 0x7fffffff) - (b[0] & 0x7fffffff)) < 0)  return -1;
		
		else
		{
			for (int i = 0; i < a_length; i++)
			{
				if      (((a[i] & 0x80000000) - (b[i] & 0x80000000)) < 0)  return  1;
				else if (((a[i] & 0x80000000) - (b[i] & 0x80000000)) > 0)  return -1;
				else if (((a[i] & 0x7fffffff) - (b[i] & 0x7fffffff)) > 0)  return  1;
				else if (((a[i] & 0x7fffffff) - (b[i] & 0x7fffffff)) < 0)  return -1;
			}
		}
	}
	
	return 0;
}


int* Math::copy(int* array, int elements)
{
	int* array1 = new int[elements];
	
	for (int i = 0; i < elements; i++)
	
	    array1[i] = array[i];
	
	return array1;
}


double Math::cos(double x)
{
	return cos_sin(x)[0];
}


double *Math::cos_sin(double x)
{

	//  The series for cosine and sine
	//
	//                n  2 n + 0
	//  cos(x) == (-1)  x       / (2n + 0)!
	//
	//                n  2 n + 1
	//  sin(x) == (-1)  x       / (2n + 1)!
	
	
	//  Divide x by 2^k
	
	int k = 0;
	
	while (x > 0.1) { x /= 2; k++; }
	
	double cosx = 0.0, sinx = 0.0;
	
	double inv_divisor = 1.0;
	
	for (int n = 0; n < 8; n++)
	{
		//  cosx += pow(x, 2*n + 0) / factorial(2*n + 0)
		//
		//      *((n % 2 == 0) ? 1 : -1);
		//
		//  (2n + 0)! == 0!, 2!, 4!, 6!, ...
		
		//  == 1, 2, 24, 720, ...
		
		if (n > 0)
		{
			inv_divisor /= 2*n - 0;
			inv_divisor /= 2*n - 1;
		}
		
		cosx += pow(x, 2*n + 0) * inv_divisor
		
		    * (((n % 2) == 0) ? 1 : -1);
	}
	
	inv_divisor = 1.0;
	
	for (int n = 0; n < 8; n++)
	{
		//  sinx += pow(x, 2*n + 1) / factorial(2*n + 1)
		//  
		//       * (((n % 2) == 0) ? 1 : -1);
		//  
		//  (2*n + 1)! == 1!, 3!, 5!, 7!, ... == 1, 6, 120, 5040, ...
		
		if (n > 0)
		{
			inv_divisor /= 2*n + 1;
			inv_divisor /= 2*n + 0;
		}
		
		sinx += pow(x, 2*n + 1) * inv_divisor
		
		    * (((n % 2) == 0) ? 1 : -1);
	}
	
	//  Multiply x by 2 ^ k using the identities
	//
	//  cos(2a) == cos^2(a) - sin^2(a),
	//
	//  sin(2a) == 2 sin(a) cos(a);
	
	for (int i = 0; i < k; i++)
	{
		double cosx1 = (cosx * cosx) - (sinx * sinx);
		
		sinx = 2 * sinx * cosx;  cosx = cosx1;
	}
	
	return new double[2] {cosx, sinx};
}


double Math::cosh(double u)
{
	//  the hyperbolic cosine function is
	//
	//                       u     -u
	//  cosh(u) = 1 / 2  ( e   +  e  )
	
	return (exp(u) + exp(-u)) / 2;
}


double *Math::cos_table(int n)
{
	double sin_table[n];
	
	double *cos_table = new double[n];
	
	int shift = n / 4;
	
	for (int i = 0; i < (n - shift); i++)
	
	    cos_table[i] = sin_table[shift + i];
	
	for (int i = 0; i < shift; i++)
	
	    cos_table[n - shift + i] = sin_table[i];
	
	return cos_table;
}


bool Math::equals(int array1[], int array2[], int elements1, int elements2)
{
	int length = elements1 <= elements2 ? elements1 :  elements2;
	
	for (int i = 0; i < length; i++)
	
	    if (array1[elements1 -1 -i] != array2[elements2 -1 -i])
	
		return false;
	
	for (int i = length; (i < elements1) && (i < elements2); i++)
	{
		if ((elements1 -1 -i) >= 0)
		
		    if (array1[elements1 -1 -i] != 0)
		
			return false;
		
		if ((elements2 -1 -i) >= 0)
		
		    if (array2[elements2 -1 -i] != 0)
		
			return false;
	}
	
	return true;
}

double Math::exp(double x)
{
	return pow(e, x);
}


vector<int> Math::expand(vector<int> vec_int, int new_size)
{
	int elements = vec_int.size();
	
	if (new_size < elements)
	{
		string message = "new size < vector size";
		
		cout << message << endl; throw string(message);
	}
	
	vector<int> vec_int1 (new_size);
	
	//  Expand the sign int
	
	int sign_int = (vec_int[0] < 0) ? -1 : 0;
	
	for (int i = 0; i < new_size; i++) vec_int1[i] = sign_int;
	
	//  Copy right to left
	
	for (int i = 0; i < elements; i++)
	
	    vec_int1[new_size -1 -i] = vec_int[elements -1 -i];
	
	return vec_int1;
}


vector<int> Math::factor(int n)
{
	//  returns a vector of prime factors
	
	int factlimit = 0x10000; // sqrt(2^32)
	
	vector<int> factors = Math::factor(n, factlimit);
	
	int factors_length = factors.size();
	
	return factors;
}


vector<int> Math::factor(int n, int max_prime)
{
	//  returns a vector of prime factors
	
	vector<int> primes = Math::primes(max_prime);
	
	int primes_length = primes.size();
	
	vector<int> factors;
	
	for (int i = 0; i < primes_length; i++)
	{
		if (is_divisible_by(n, primes[i]))
		{
			factors.push_back(int(primes[i]));
			
			n /= primes[i--];
		}
	}
	
	return factors;
}


vector<int> Math::factor(Number& number, int max_prime)
{
	//  returns a vector of prime factors less than prime
	
	Number n (number);
	
	if (!n.is_integer()) throw string(
	
	    "factorization is only defined for integers");
	
	vector<int> *list = new vector<int>();
	
	vector<int> primes = Math::primes(max_prime);
	
	int primes_length = primes.size();
	
	if (n.signum() == -1) list->push_back(-1);
	
	for (int i = 0; i < primes_length; i++)
	{
		if (n.is_divisible_by(primes[i]))
		{
		      list->push_back(primes[i]);
		    
		      n = n.divide(primes[i--]);
		}
	}
	
	return *list;
}



void Math::flip_bit(vector<int> a, long bit)
{
	int a_length = a.size();
	
	int int1 = a[a_length -1 - (int) (bit / 32)];
	
	int1 ^= (1 << (bit % 32));
	
	a[a_length -1 - (int) (bit / 32)] = int1;
}

void Math::flip_bit(int array[], int elements, long bit)
{
	int array_length = elements;
	
	int int1 = array[array_length -1 - (int) (bit / 32)];
	
	int1 ^= (1 << (bit % 32));
	
	array[array_length -1 - (int) (bit / 32)] = int1;
}


long Math::gcd(long m, long n)
{
	if ((m == 0) || (n == 0))
	
	    return 1;
	
	m = abs(m);  n = abs(n);
	
	long a, b, c = 1;
	
	if (m >= n)
	
	     { a = m; b = n; }
	else { a = n; b = m; }
	
	do { c = a % b; a = b; b = c; }
	
	while (c != 0);
	
	return a;
}


int Math::get_bit(int array[], int elements, long bit)
{
	if ((bit / 32) >= elements)
	{
		string message = "bit / 32 > array size";
		
		cout << message << endl; throw string(message);
	}
	
	int int_32 = array[elements -1 - (int) (bit / 32)];
	
	unsigned int bit_mask[32] =
	{
	      0x00000001, 0x00000002, 0x00000004, 0x00000008,
	      0x00000010, 0x00000020, 0x00000040, 0x00000080,
	      0x00000100, 0x00000200, 0x00000400, 0x00000800,
	      0x00001000, 0x00002000, 0x00004000, 0x00008000,
	      0x00010000, 0x00020000, 0x00040000, 0x00080000,
	      0x00100000, 0x00200000, 0x00400000, 0x00800000,
	      0x01000000, 0x02000000, 0x04000000, 0x08000000,
	      0x10000000, 0x20000000, 0x40000000, 0x80000000,
	};
	
	return  (int_32 & bit_mask[bit % 32]) == 0 ? 0 : 1;
}

int Math::get_bit(vector<int>& vec_int, long bit)
{
	int elements = vec_int.size();
	
	if ((bit / 32) >= elements)
	{
		string message = "bit / 32 > vector size";
		
		cout << message << endl; throw string(message);
	}
	
	int int_32 = vec_int[elements -1 - (int) (bit / 32)];
	
	unsigned int bit_mask[32] =
	{
	      0x00000001, 0x00000002, 0x00000004, 0x00000008,
	      0x00000010, 0x00000020, 0x00000040, 0x00000080,
	      0x00000100, 0x00000200, 0x00000400, 0x00000800,
	      0x00001000, 0x00002000, 0x00004000, 0x00008000,
	      0x00010000, 0x00020000, 0x00040000, 0x00080000,
	      0x00100000, 0x00200000, 0x00400000, 0x00800000,
	      0x01000000, 0x02000000, 0x04000000, 0x08000000,
	      0x10000000, 0x20000000, 0x40000000, 0x80000000,
	};
	
	return  (int_32 & bit_mask[bit % 32]) == 0 ? 0 : 1;
}


int Math::get_lowest_set_bit(long a)
{
	if (a == 0) return -1;
	
	int bit = 0;
	
	while ((a & 1) == 0)
	
	    { a >>= 1; bit++; }
	
	return bit;
}


long Math::get_lowest_set_bit(vector<int>& a)
{
	long bit = 0;
	
	int array_length = a.size();
	
	while ((test_bit(a, bit) == false)
	
	    && (bit < 32 * array_length))  bit++;
	
	if (bit >= 32 * array_length) return -1;
	
	return bit;
}


long Math::get_lowest_set_bit(int array[], int elements)
{
	long bit = 0;
	
	int array_length = elements;
	
	while ((test_bit(array, array_length, bit) == false)
	
	    && (bit < 32 * array_length))  bit++;
	
	if (bit >= 32 * array_length) return -1;
	
	return bit;
}


double Math::inverse(double n)
{
	return 1 / n;
}

bool Math::is_coprime_with(long m, long n)
{
	return gcd(m, n) == 1;
}

bool Math::is_divisible_by(long m, long n)
{
	return (m % n) == 0 ? true : false;
}

bool Math::is_power(long n, int exp)
{
	int root = (int) (Math::root(n, exp) + 0.1);
	
	return (pow(root, exp) - n) == 0.0 ? true : false;
}

bool Math::is_power_of_2(long n)
{
	return is_power_of(2, n);
}

bool Math::is_power_of(int base, long n)
{
	int b = base;
	
	if (b == 2) return
	
	    ((n & (n - 1)) == 0) ? true : false;
	
	long q = n / base;
	
	int qsize = bit_count(q);
	int bsize = bit_count(b);
	
	int dsize = qsize / bsize;
	
	long d = (long) pow(b, dsize-2);
	
	if ((n % d) != 0) return false;
	
	q /= d;
	
	while (q != 1)
	{
		q /= b;
		
		if (q == b)
		
		    return true;
	}
	
	return false;
}



bool Math::is_prime(int n)
{

	int primes[168] =
	{
	       2,  3,    5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,
	      59,  61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131,
	     137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
	     227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
	     313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
	     419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
	     509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613,
	     617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719,
	     727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827,
	     829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
	     947, 953, 967, 971, 977, 983, 991, 997
	};
	
	
	//  Divide by small primes
	
	int primes_size = sizeof(primes) / sizeof(primes[0]);
	
	for (int i = 0; i < primes_size; i++)
	{
		if ((n % primes[i]) == 0)
		{
			if (n > primes[primes_size -1])
			
			    return false;
			
			else if (n / primes[i] == 1)
			
			    return true;
		}
		
		else  return false;
	}
	
	
	//  Apply the Fermat test a^(n-1) mod n == 1
	//  using bases a = 2 and a = 3 as witnesses
	
	//  If n is prime, then the totient is phi(n) == n-1.
	
	if ( (mod_pow(2, n-1, n) != 1)
	  || (mod_pow(3, n-1, n) != 1) )
	
	     return false; // composite
	
	
	//  Test for pseudoprimes using the Miller-Rabin test
	
	//  ...   ...
	
	//  ...   ...
	
	//  ...   ...
	
	
	return true;
}


bool Math::is_quadratic_residue(int residue, int p)
{	
	if (!is_prime(p))
	{
		string message = "non-prime modulus";
		
		cout << message << endl; throw string(message);
	}
	
	if (residue == 0) return false;
	
	return mod_pow(residue, phi(p) / 2, p) == 1 ? 
	
	    true : false;
}


bool Math::is_square(double val)
{
	if (val == 0.0) return true;
	
	if (square((long) sqrt(val)) == val)
	
	     return true;
	
	else return false;
}


bool Math::is_square(long val)
{
	return (square((long) sqrt(val)) == val) ?
	
	    true : false;
}


int Math::lambda(int n)
{

	//  Carmichael lambda function
	//
	//  a^lambda(n) == 1 mod(n) for all (a, n) == 1
	//
	//  lambda(n) == lcm(phi(fact(n)))
	
	
	//  Factor the modulus
	
	vector<int> factors = factor(n);
	
	int factors_length = factors[0];
	
	//  Calculate the lcm of the phi functions recursively
	
	int lcm = 1;
	
	for (int i = 1; i <= factors_length; i++)
	
	    lcm = (int) Math::lcm(lcm, phi(factors[i]));
	
	return lcm;
}


long Math::lcm(int a, int b)
{
	//  computes the least common multiple
	
	//  a x b == gcd(a, b) x lcm(a, b)
	
	return 1L * a * b / gcd(a, b);
}


Number& Math::lcr(int *r, int *n, int elements)
{
	int size = elements;
	
	Number r1[size], n1[size];
	
	for (int i = 0; i < size; i++)
	{
		r1[i] = Number(r[i]).lvalue();
		n1[i] = Number(n[i]).lvalue();
	}
	
	return lcr(r1, n1, size);
}


Number& Math::lcr(Number r[], Number n[], int elements)
{

	//  computes the least common / composite / Chinese
	//  remainder or residue from a set of reduced
	//  residues and moduli
	
	
	int r_length = elements;
	int n_length = elements;
	
	
	//  t = the number of coprime moduli
	
	int t = n_length;
	
	
	//  For each n[i], compute the inverse of n[i] modulo p[i]
	
	Number m[t];
	
	for (int i = 0; i < t; i++)
	{
		Number p = n[i];
		
		Number N(1);
		
		for (int j = 0; j < t; j++)
		
		    if (j != i)  N = N.multiply(n[j]);
		
		N = N.mod(p);
		
		if (!N.equals(0))
		
		    m[i] = N.mod_inverse(p);
		
		else // p[i] == 1, m[i] = any number, r[i] = 0
		
		    m[i] = Number(0).lvalue(); // r[i] mod 1 equals 0
	}
	
	
	//  Compute Ni, Mi, and the composite residue
	//
	//  R = the sum over i of r[i] N[i] M[i] (mod N)
	
	
	Number R (0);
	
	
	for (int i = 0; i < t; i++)
	{
		Number N (1);
		
		for (int j = 0; j < t; j++)
		
		    if (j != i)  N = N.multiply(n[j]);
		
		R = R .add( r[i] .multiply(m[i]) . multiply(N) );
	}
	
	
	//  Compute the composite modulus N
	
	Number N(1);
	
	for (int i = 0; i < t; i++)
	
	    N = N .multiply(n[i]);
	
	
	//  Reduce the common residue to its least value modulo N
	
	return  R .mod(N);
}


double Math::log(double a, double b)
{
	//  returns log(a) to the base b
	
	return log(a) / log(b);
}


double Math::log(double d)
{

	//  computes the natural logarithm
	//
	//  ln(x) for any real x > 0
	
	
	//  The series for log(1 + x) is
	//
	//                   n+1   n
	//  log(1 + x) = (-1)    x  /  n
	//
	//  for  -1 < x < 1  or  0 < (1 + x) < 2
	
	
	//  Set d = e^k * b. Then
	//
	//  log(d) == log(e^k * b) == k + log(b)
	
	int k = 0;
	
	double d1 = d;
	
	while (d1 > 1.5) { d1 /= e; k++; }
	
	double b = d1, x = b - 1, log = k;
	
	for (long n = 1; n < 32; n++)
	{
		double t = pow(x, n) / n;
		
		if ((n % 2) == 1)
		
		     log += t;
		
		else log -= t;
	}
	
	//  Verify that  e^log(x) / x == x / x == 1
	
	double q = abs(exp(log) / d);
	
	if ((q > 1.00001) || (q < 0.9999))
	{
		cout << string("d == ") << d
		
		    << string(" log(d) == ") << log << endl;
		
		cout << string("e^log(d) / d == ")
		
		    << (exp(log) / d) << endl;
		
		throw string("arithmetic exception");
	}
	
	return log;
}


double Math::log(int a, int b)
{
	return log((double) a, (double) b);
}

double Math::log(int a)
{
	return log((double) a);
}

float Math::log(float a, float b)
{
	return (float) log((double) a, (double) b);
}

float Math::log(float f)
{
	return (float) log((double) f);
}

int Math::log2(int n)
{
	//  returns the bit count-1
	
	return bit_count(n) - 1;
}

double Math::log2(double n)
{
	//  computes the log to the base 2
	//
	//  using the logarithmic identity
	
	//  log n == log n / log b
	//     b        x       x
	//
	//  where x can be any base
	
	
	return log(n) / log(2);
}


double Math::log10(double n)
{
	//  computes the common logarithm
	
	return log(n) / log(10);
}


int Math::mod(vector<int> vec_int, int val)
{
	//  reduces a vector modulo an int
	
	//  This method is much faster for computing residues
	//  than using quadratic division because it computes
	//  the residue in O(n) steps or operations.
	
	//  Compute the mods of the ints and add them
	
	long residue = 0L;
	
	long twosr = 1L; // twos residue
	
	for (int i = vec_int.size() -1; i >= 0; i--)
	{
		long a = vec_int[i] & 0xffffffffL;
		
		long r = a % val;
		
		r = (r * twosr) % val;
		
		twosr *= 0x100000000L;
		
		twosr %= val;
		
		residue += r;
	}
	
	residue %= val;
	
	//  This test was used for debugging
	//
	//  if (!Number(vec_int) .mod(Number(val).lvalue()) .equals(residue))
	//  {
	//	string message = "mod(vector, int) error";
	//	
	//	cout << message << endl; throw string(message);
	//  }
	
	return (int) residue;
}


int Math::mod_inverse(int a1, int n)
{
	//  computes the mod inverse using the extended Euclidean algorithm
	
	//  The int mod_inverse method is much faster
	//  than the Number mod_inverse method
	
	//  variables:  a, b, c, q, r,  x, x1, x2,  y, y1, y2;
	
	if (a1 < 0)  a1 = ((a1 % n) + n) % n;
	
	long a, b, c, q, r;
	
	a = b = c = q = r = 0;
	
	long x, x1, x2, y, y1, y2;
	
	x1 = 0;  x2 = 1;
	y1 = 1;  y2 = 0;
	
	a = abs(a1 % n); b = n;
	
	while (a != 0)
	{
		q = b / a;
		
		r = b - (q * a);
		
		x = x2 - (q * x1);
		y = y2 - (q * y1);
		
		x2 = x1;  x1 = x;
		y2 = y1;  y1 = y;
		
		b = a;  a = r;
	}
	
	c = b; x = x2; y = y2;
	
	long inva = y;
	
	if (c != 1)
	{
		string message = "non-invertible number";
		
		cout << message << endl; throw string(message);
	}
	
	//  Check the answer
	
	long product = (((inva * (a1 % n)) % n) + n) % n;
	
	if (product != 1)
	{
		string message = "modular inversion error";
		
		cout << message << endl; throw string(message);
	}
	
	inva %= n;
	
	if (inva < 0) inva += n;
	
	return (int) (inva % n);
}



int Math::mod_pow(int base, int exp, int mod)
{
	int n = mod,  x = exp;
	
	long a = base % n;
	
	long y = 1L;
	
	bool modis_power_of_2;
	
	if (((mod - 1) & mod) == 0)
	
	     modis_power_of_2 = true;
	else modis_power_of_2 = false;
	
	int modm1 = mod - 1;
	
	//  Use the square and multiply method
	
	while (x != 0)
	{
		//  Accumulate squares
		
		if ((x % 2) == 1)
		{
			y *= a;
			
			if (modis_power_of_2)
			
			     y &= modm1;
			
			else y %= n;
		}
		
		//  Square the square
		
		a *= a;
		
		if (modis_power_of_2)
		
		     a &= modm1;
		
		else a %= n;
		
		//  Shift the exponent to the next bit
		
		x >>= 1;
	}
	
	if ((base < 0) && (exp % 2 == 0))
	
	    y = abs(y);
	
	
	y = ((y % mod) + mod) % mod;
	
	return (int) y;
}



vector<int> Math::multiply(vector<int> multiplier, vector<int> multiplicand)
{

	//  multiplies two vectors
	
	//  Count the number of twos or zeros
	
	int twos1 = 0, twos2 = 0;
	
	int multiplier_length   = multiplier.size();
	int multiplicand_length = multiplicand.size();
	
	for (int i = 0; i < multiplier_length; i++)
	{
		if (multiplier[multiplier_length -1 -i] == 0) twos1++;
		
		else break;
	}
	
	for (int i = 0; i < multiplicand_length; i++)
	{
		if (multiplicand[multiplicand_length -1 -i] == 0) twos2++;
		
		else break;
	}
	
	
	//  Right shift the multiplier and multiplicand to remove the trailing zeros
	//  and then trim to remove the leading zeros
	
	vector<int> multiplier1   = shift_right(multiplier,   32 * twos1);
	vector<int> multiplicand1 = shift_right(multiplicand, 32 * twos2);
	
	multiplier1   = trim(multiplier1,   twos1);
	multiplicand1 = trim(multiplicand1, twos2);
	
	
	//  Compute the product using the quadratic multiplier
	
	vector<int> product1 = quad_multiply(multiplier1, multiplicand1);
	
	
	//  Replace the twos or zeros in the product
	
	int twos = twos1 + twos2;
	
	product1 = shift_left(product1, 32*twos, 32*twos);
	
	
	//  Make a new copy of the vector
	
	vector<int> product = product1;
	
	return product;
}


vector<int> Math::multiply(vector<int> vec_int, int y)
{
	int m1 = vec_int.size(), x[m1];
	
	for (int i = 0; i < m1; i++) x[i] = 0;
	for (int i = 0; i < m1; i++) x[i] = vec_int[i];
	
	int *w = new int[m1 + 1];
	
	for (int i = 0; i < (m1 + 1); i++) w[i] = 0;
	
	int z = 0;
	
	unsigned long uv = 0L;
	
	long xy = 0L;
	
	for (int j = 0; j < m1; j++)
	{
		//  Compute the product of x[j] and y
		
		xy = 1L * x[m1 -1 -j] * y;
		
		//  Convert the signed product to an unsigned product
		
		if (x[m1 -1 -j] < 0) xy += ((1L * y)           << 32);
		if           (y < 0) xy += ((1L * x[m1 -1 -j]) << 32);
		
		//  //  Add the product and carry bits to the sum
		
		uv = xy + (w[(m1 + 1) -1 -j] & 0xffffffffL) + (z & 0xffffffffL);
		
		//  Copy the lower int of uv to the product array
		
		w[(m1 + 1) -1 -j] = uv;
		
		//  Copy the upper int of uv to the carry int z
		
		z = ((uv >> 32) & 0xffffffffL);
	}
	
	//  Copy the upper int of uv to w
	
	w[(m1 + 1) -1 -m1] = ((uv >> 32) & 0xffffffffL);
	
	vector<int> product (m1 + 1);
	
	for (int i = 0; i < m1 + 1; i++)
	
	    product[i] = w[i];
	
	delete[] w;
	
	return product;
}


int Math::phi(int n)
{

	if ((n == 1) || (n == 2)) return 1;
	
	if (n == 3) return 2;
	
	if (is_prime(n)) return n - 1;
	
	
	int factlimit = 0x10000;
	
	vector<int> factors = Math::factor(n, factlimit);
	
	int size = factors.size();
	
	//  Sort and collate the prime factors
	
	//  (First convert the vector to an array to use the int[] sorter)
	
	int array[size]; for (int i = 0; i < size; i++) array[i] = factors[i];
	
	vector<int*> int_pow = Math::sort_and_collate(array, size);
	
	int int_pow_length = int_pow[0][0];
	
	//  Compute the totient using the list of prime powers
	
	int totient = 1;
	
	for (int i = 0; i < int_pow_length; i++)
	
	    totient = totient * (int_pow[i][0] - 1) *
	
		(int) pow(int_pow[i][0], int_pow[i][1] - 1);
	
	
	//  Verify the totient using the formula
	//
	//  a ^ phi(n) == 1 (mod n)  where (a, n) == 1
	
	//  Choose a random base a
	
	int a = n / 2;
	
	if (a < 2) a = 2;
	
	while (gcd(a, n) != 1) a++;
	
	//  if (mod_pow(a, totient, n) != 1)
	//
	//	throw exception
	
	return totient;
}


double Math::pow(double base, int exp)
{
	return pow(base, (long) exp);
}

double Math::pow(double base, long exp)
{
	if (exp == 1) return base;
	
	double a = base, y = 1.0;
	
	long b = exp,  x = abs(b);
	
	//  Use the square and multiply method
	
	while (x != 0)
	{
		//  Accumulate squares
		
		if ((x & 1) == 1)
		
		    y *= a;
		
		//  Square the square
		
		a *= a;
		
		//  Shift the exponent to the next bit
		
		x >>= 1;
	}
	
	//  Invert the output if the exponent is negative
	
	if (b < 0)  y = 1 / y;
	
	return y;
}


double Math::pow(double base, double exp)
{

	//  computes the function
	//
	//          x     log a  x      k x
	//  y  =  a  == (e     )   == e
	//  
	//  for a real base and exponent
	//  
	//  (integer exponents should use pow(double, long))
	
	
	double a = base, b = exp;
	
	if (b == 0) return 1.0;
	
	
	//  If the exponent is an integer return y_int
	
	if (abs(b - (int) (b)) < 0.00000001)
	
	    return pow(a, (long) b);
	
	
	//  Solve for k in the equation
	//         k
	//  a == e  and multiply b (= exp) by k
	
	double k = (abs(a - e) > 0.00000001)? k = log(a) : 1;
	
	//            x     k x    k x   x'
	//  Compute a  == (e ) == e == e
	//
	//  where k = log(a) from the series
	//
	//   x      n
	//  e  == x  / n!
	//
	//  where x is the exponent and
	//  n is the index 0, 1, 2, ...
	//
	//  This series brings the fraction x from
	//  the exponent down to the base to avoid
	//  doing square roots which are expensive.
	
	
	double y_frac = 0.0, temp = 1.0;
	
	//  temp = x^0 / 0!
	
	double factorial = 1.0;
	
	//  Compute the new exponent x' = k x
	
	double b1 = b * k;
	
	//  Compute the integer value y_int
	
	double y_int = pow(e, (long) b1);
	
	//  Compute the fractional value y_frac
	
	b1 = abs(b1 - (long) b1);
	
	for (int n = 1; n < 32; n++)
	{
		y_frac += (temp / factorial);
		
		temp *= b1;
		
		factorial *= n;
	}
	
	//  Compute y = y_int * y_frac
	
	double y = y_int * y_frac;
	
	//  Invert the output if the exponent is negative
	
	if (b <= 0)
	
	    y = inverse(y);
	
	return y;
}


vector<int> Math::primes(int n)
{

	//  returns the first n primes using the prime number sieve
	
	//  The prime number sieve works by striking out integers that
	//  are multiples of integers starting from the number two.
	//
	//  The prime number sieve is thousands of times faster than the
	//  Fermat primality test for finding the first few million primes.
	
	//  Example  vector<int> primes = Math::primes(168);
	//
	//  [ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
	//   53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
	//   113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
	//   181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241,
	//   251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
	//   317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389,
	//   397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
	//   463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
	//   557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617,
	//   619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691,
	//   701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773,
	//   787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859,
	//   863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
	//   953, 967, 971, 977, 983, 991, 997 ]
	
	
	int size = n + 1;
	
	bool sieve[size];
	
	sieve[0] = false;
	sieve[1] = false;
	
	for (int i = 2; i < size; i++)
	
	    sieve[i] = true;
	
	for (int i = 0; i < size; i++)
	{
		if (sieve[i] == false) continue;
		
		for (int j = 2; i*j < size; j++)
		
		    sieve[i*j] = false;
	}
	
	vector<int> list;
	
	for (int i = 0; i < size; i++)
	
	    if (sieve[i] == true)
	
		list.push_back(i);
	
	return list;
}



vector<int> Math::quad_multiply(vector<int> vec_int1, vector<int> vec_int2)
{
	//  multiplies two vectors using a quad multiplier
	//
	//  The length of the product vector equals the sum of the
	//  lengths of the two vectors.
	
	
	int m1 = vec_int1.size();
	int m2 = vec_int2.size();
	
	//  If the multiplier and multiplicand are the same,
	//  use the quadratic squarer instead of the multiplier
	//  because squaring is twice as fast as multiplying
	
	if (m1 == m2)
	{
		bool equals = true;
		
		for (int i = 0; i < m1; i++)
		
		    if (vec_int1[i] != vec_int2[i])
		
			{ equals = false; break; }
		
		////  if (equals) return quad_square(x);
	}
	
	
	int x[m1], y[m2];
	
	for (int i = 0; i < m1; i++) x[i] = 0;
	for (int i = 0; i < m2; i++) y[i] = 0;
	
	for (int i = 0; i < m1; i++) { x[i] = vec_int1[i]; }
	for (int i = 0; i < m2; i++) { y[i] = vec_int2[i]; }
	
	int w[m1 + m2];
	
	for (int i = 0; i < (m1 + m2); i++) w[i] = 0;
	
	for (int i = 0; i < m2; i++)
	{
		int  z = 0;
		
		unsigned long uv = 0L;
		
		long xy = 0L;
		
		for (int j = 0; j < m1; j++)
		{
			//  Compute the unsigned product of x[j] and y[i]
			
			xy = 1L * x[m1 -1 -j] * y[m2 -1 -i];
			
			//  Convert the signed product to an unsigned product
			
			   if (x[m1 -1 -j] < 0) xy += ((1L * y[m2 -1 -i]) << 32);
			   if (y[m2 -1 -i] < 0) xy += ((1L * x[m1 -1 -j]) << 32);
			
			//  Add the product and carry bits to the sum
			
			uv = xy + (w[(m1 + m2) -1 -(i + j)] & 0xffffffffL) + (z & 0xffffffffL);
			
			//  Copy the lower int of uv to the product array w
			
			w[(m1 + m2) -1 - (i + j)] = uv;
			
			//  Copy the upper int of uv to the carry int z
			
			z = ((uv >> 32) & 0xffffffffL);
		}
		
		//  Copy the upper int of uv to w
		
		w[(m1 + m2) -1 - (i + m1)] = ((uv >> 32) & 0xffffffffL);
	}
	
	vector<int> product (m1 + m2);
	
	for (int i = 0; i < m1 + m2; i++)
	
	    product[i] = w[i];
	
	return product;
}


int *Math::quad_multiply(int array1[], int array2[], int elements1, int elements2)
{
	//  multiplies two arrays using a quad multiplier
	//
	//  The length of the product array equals the sum of the
	//  lengths of the two arrays.
	
	int m1 = elements1;
	int m2 = elements2;
	
	//  If the multiplier and multiplicand are the same,
	//  use the quadratic squarer instead of the multiplier
	//  because squaring is twice as fast as multiplying
	
	if (m1 == m2)
	{
		bool equals = true;
		
		for (int i = 0; i < m1; i++)
		
		    if (array1[i] != array2[i])
		
			{ equals = false; break; }
		
		////  if (equals) return quad_square(x);
	}
	
	
	int x[m1], y[m2];
	
	for (int i = 0; i < m1; i++) x[i] = 0;
	for (int i = 0; i < m2; i++) y[i] = 0;
	
	for (int i = 0; i < m1; i++) { x[i] = array1[i]; }
	for (int i = 0; i < m2; i++) { y[i] = array2[i]; }
	
	int *w = new int[m1 + m2];
	
	for (int i = 0; i < (m1 + m2); i++) w[i] = 0;
	
	for (int i = 0; i < m2; i++)
	{
		int  z = 0;
		
		unsigned long uv = 0L;
		
		long xy = 0L;
		
		for (int j = 0; j < m1; j++)
		{
			//  Compute the unsigned product of x[j] and y[i]
			
			xy = 1L * x[m1 -1 -j] * y[m2 -1 -i];
			
			//  Convert the signed product to an unsigned product
			
			   if (x[m1 -1 -j] < 0) xy += ((1L * y[m2 -1 -i]) << 32);
			   if (y[m2 -1 -i] < 0) xy += ((1L * x[m1 -1 -j]) << 32);
			
			//  Add the product and carry bits to the sum
			
			uv = xy + (w[(m1 + m2) -1 -(i + j)] & 0xffffffffL) + (z & 0xffffffffL);
			
			//  Copy the lower int of uv to the product array w
			
			w[(m1 + m2) -1 - (i + j)] = uv;
			
			//  Copy the upper int of uv to the carry int z
			
			z = ((uv >> 32) & 0xffffffffL);
		}
		
		//  Copy the upper int of uv to w
		
		w[(m1 + m2) -1 - (i + m1)] = ((uv >> 32) & 0xffffffffL);
	}
	
	return w;
}


double Math::root(double n, int k)
{
	//  computes the kth root of a number n
	
	if (n == 0) return 0;
	if (k == 1) return n;
	
	if ((k < 1) || ((k % 2 == 0) && (n < 0)))
	
	    throw string("illegal argument exception");
	
	if (n == 0)  return 0;
	if (k == 1)  return n;
	
	//  Set r approximately equal to root(n)
	
	double r = 1.0;
	
	if (r < n) while (true)
	{
		double pow = 1.0;
		
		for (int i = 0; i < k; i++) pow *= r;
		
		if (pow < n)  r *= 16*16*16;
		
		else break;
	}
	
	else if (r > n) while (true)
	{
		double pow = 1.0;
		
		for (int i = 0; i < k; i++) pow *= r;
		
		if (pow > n)  r /= 16*16*16;
		
		else break;
	}
	
	for (int j = 0; j < 16; j++)
	{
		double pow = Math::pow(r, k-1);
		
		r = ((k-1) * r  +  n / pow) / k;
	}
	
	if (abs(pow(r, k) / n) > 1.0001)
	{
		cout << "radicand == " << n;
		cout << "    root == " << r;
		cout << "| root^exp | == " << abs(pow(r, k));
		cout << "root^exp - radicand == " << abs(pow(r, k) - n);
		
		throw string("root extraction error");
	}
	
	return r;
}


void Math::set_bit(vector<int>& a, long bit)
{
	a[a.size() - 1 - (int) (bit / 32)] |= (1 << (bit % 32));
}

void Math::set_bit(int array[], int elements, long bit)
{
	array[elements - 1 - (int) (bit / 32)] |= (1 << (bit % 32));
}


vector<int> Math::shift_left(vector<int> vec_int, long bits)
{
	int vec_length = vec_int.size();
	
	if ((bits > 32 * vec_length) || (bits < 0))
	{
		string message = "0 <= bits < vector size";
		
		cout << message << endl; throw string(message);
	}
	
	int ints = (int) (bits / 32);
	int b    = (int) (bits % 32);
	
	vector<int> vec_int1 = vec_int;
	
	if (ints > 0)
	{
		for (int i = 0; (i + ints) < vec_length; i++)
		
		    vec_int1[i] = vec_int[i + ints];
		
		for (int i = vec_length - ints; i < vec_length; i++)
		
		    vec_int1[i] = 0;
	}
	
	long carry = 0L;
	
	for (int i = vec_length -1; i >= 0; i--)
	{
		long a = (vec_int1[i] & 0xffffffffL);
		
		vec_int1[i] = (int) (((a << b) + carry) & 0xffffffffL);
		
		carry = a >> (32 - (1L * b));
	}
	
	return vec_int1;
}


int *Math::shift_left(int array[], int elements, long bits)
{
	int array_length = elements;
	
	if ((bits > 32 * array_length) || (bits < 0))
	{
		string message = "0 <= bits < array size";
		
		cout << message << endl; throw string(message);
	}
	
	int ints = (int) (bits / 32);
	int b    = (int) (bits % 32);
	
	int *array1 = new int[array_length];
	
	for (int i = 0; i < array_length; i++)
	
	    array1[i] = array[i];
	
	if (ints > 0)
	{
		for (int i = 0; (i + ints) < array_length; i++)
		
		    array1[i] = array[i + ints];
		
		for (int i = array_length - ints; i < array_length; i++)
		
		    array1[i] = 0;
	}
	
	long carry = 0L;
	
	for (int i = array_length -1; i >= 0; i--)
	{
		long a = (array1[i] & 0xffffffffL);
		
		array1[i] = (int) (((a << b) + carry) & 0xffffffffL);
		
		carry = a >> (32 - (1L * b));
	}
	
	return array1;
}


vector<int> Math::shift_left(vector<int> vec_int, long expansion, long bits)
{
	//  expands the vector before left shifting
	
	if (bits < 0)
	{
		string message = "bits < 0";
		
		cout << message << endl; throw string(message);
	}
	
	int vec_length = vec_int.size();
	
	int vec_length1 = vec_length + ((expansion + 31) / 32);
	
	vector<int> vec_int1 (vec_length1);
	
	for (int i = 0; i < vec_length1; i++) vec_int1[i] = 0;
	
	//  Copy right (to left)
	
	for (int i = 0; i < vec_length; i++)
	
	     vec_int1[vec_length1 -1 -i]
	   = vec_int [vec_length  -1 -i];
	
	vector<int> vec_int2 = shift_left(vec_int1, bits);
	
	return vec_int2;
}


int *Math::shift_left(int array[], int elements, long expansion, long bits)
{
	//  expands the array before left shifting
	
	if (bits < 0)
	{
		string message = "bits < 0";
		
		cout << message << endl; throw string(message);
	}
	
	int array_length = elements;
	
	int array_length1 = array_length + ((expansion + 31) / 32);
	
	int *array1 = new int[array_length1];
	
	for (int i = 0; i < array_length1; i++) array1[i] = 0;
	
	//  Copy right (to left)
	
	for (int i = 0; i < array_length; i++)
	
	     array1[array_length1 -1 -i]
	   = array [array_length  -1 -i];
	
	int *array2 = shift_left(array1, array_length1, bits);
	
	delete[] array1;
	
	return array2;
}


vector<int> Math::shift_right(vector<int> vec_int, long bits)
{
	if (bits < 0)
	{
		string message = "bits < 0";
		
		cout << message << endl; throw string(message);
	}
	
	int vec_length = vec_int.size();
	
	if (bits > 32 * vec_length) bits = 32 * vec_length;
	
	int ints = (int) (bits / 32);
	int b    = (int) (bits % 32);
	
	vector<int> vec_int1 (vec_length);
	
	for (int i = 0; i < vec_length; i++)
	
	    vec_int1[i] = vec_int[i];
	
	if (ints > 0)
	{
		for (int i = 0; (i + ints) < vec_length; i++)
		
		      vec_int1[vec_length -1 - i]
		    = vec_int [vec_length -1 - (i + ints)];
		
		for (int i = 0; i < ints; i++)
		
		    vec_int1[i] = 0;
	}
	
	long carry = 0L;
	
	for (int i = 0; i < vec_length; i++)
	{
		unsigned long a = (vec_int1[i] & 0xffffffffL);
		
		vec_int1[i] = (((a >> b) + carry) & 0xffffffffL);
		
		carry = ( a << ( 32 - (1L * b)) );
	}
	
	return vec_int1;
}


int *Math::shift_right(int array[], int elements, long bits)
{
	int array_length = elements;
	
	if ((bits > 32 * array_length) || (bits < 0))
	{
		string message = "0 <= bits < array size";
		
		cout << message << endl; throw string(message);
	}
	
	int ints = (int) (bits / 32);
	int b    = (int) (bits % 32);
	
	int *array1 = new int[array_length];
	
	for (int i = 0; i < array_length; i++)
	
	    array1[i] = array[i];
	
	if (ints > 0)
	{
		for (int i = 0; (i + ints) < array_length; i++)
		
		      array1[array_length -1 - i] =
		      array [array_length -1 - (i + ints)];
		
		for (int i = 0; i < ints; i++)
		
		    array1[i] = 0;
	}
	
	long carry = 0L;
	
	for (int i = 0; i < array_length; i++)
	{
		unsigned long a = (array1[i] & 0xffffffffL);
		
		array1[i] = (((a >> b) + carry) & 0xffffffffL);
		
		carry = ( a << ( 32 - (1L * b)) );
	}
	
	return array1;
}


double Math::sin(double x)
{
	return cos_sin(x)[1];
}

double *Math::sin_table(int n)
{
	if (n == 0) return nullptr;
	
	double *sin = new double[n];
	
	double cos[n];
	
	if ((n % 2) == 0)
	{
		sin[0] = 0.0;  cos[0] = 1.0;
		
		if ((n == 1) || (n == 2)) return sin;
		
		sin[1] = Math::sin(2 * 3.141592653589793 / n);
		
		for (int i = 2; i <= n / 4; i++)
		
		    sin[i] = Math::sin(2 * 3.14159265358979 * i / n);
		
		for (int i = 0; i <= n / 4; i++)
		
		    sin[1*n / 4 + i] = sin[1*n / 4 - i];
		
		for (int i = 0; i < n / 2; i++)
		
		    sin[2*n / 4 + i] = -sin[0*n / 4 + i];
	}
	
	else for (int i = 0; i < n; i++)
	
	    sin[i] = Math::sin(2 * 3.14159265358979 * i / n);
	
	return sin;
}


Number *Math::sin_table(Number& n)
{
	int m = n.int_value();
	
	Number *sin = new Number[m];
	
	if ((m % 2) == 0)
	{
		sin[0] = Number(0).lvalue();
		
		if (m == 1) return sin;
		
		int precision = n.get_precision() == 0 ? 8 : n.get_precision();
		
		sin[1] = Number::pi(precision).multiply(2).divide(n).sin();
		
		if (m == 2) return sin;
		
		for (int i = 2; i <= m / 4; i++)
		
		    sin[i] = Number::pi(precision).multiply(2)
		
			.divide(n).multiply(i).sin();
		
		for (int i = 0; i <= m / 4; i++)
		
		    sin[1*m / 4 + i] = sin[1*m / 4 - i];
		
		for (int i = 0; i < m / 2; i++)
		
		    sin[2*m / 4 + i] = sin[0*m / 4 + i] .negate();
	}
	
	else
	{	int precision = n.get_precision() == 0 ? 8 : n.get_precision();
		
		for (int i = 0; i < m; i++)
		
		    sin[i] = Number::pi(precision).multiply(2)
		
			.divide(n).multiply(i).sin();
	}
	
	return sin;
}


double Math::sinh(double u)
{
	//  The hyperbolic sine function is
	//
	//                        u     -u
	//  sinh(x)  =  1 / 2  ( e  -  e  )
	
	return (exp(u) - exp(-u)) / 2;
}



void Math::sort(vector<int> & a)
{
	int elements = a.size();
	
	int array[elements];
	
	for (int i = 0; i < elements; i++)
	
	    array[i] = a[i];
	
	sort(array, elements);
	
	vector<int> b;
	
	for (int i = 0; i < elements; i++)
	
	    b.push_back(array[i]);
	
	a = b;
}


void Math::sort(vector<vector<int>>& a)
{
	//  Copy and pad the vector
	
	vector<vector<int>> x = a;
	
	//  Find the smallest element
	
	int s = 0;
	
	for (int i = 0; i < a.size(); i++)
	{
		int int1 = a[i][0];
		
		if (int1 < s) s = int1;
	}
	
	if ((x.size() & (x.size() -1)) != 0)
	{
		//  Pad the vector to a power of 2
		
		int size = x.size();
		
		int power_of_two = 1;
		
		while (size != 0)
		
		    { size >>= 1; power_of_two <<= 1; }
		
		vector<vector<int>> vec_int (power_of_two);
		
		for (int i = 0; i < x.size(); i++)
		
		    vec_int[vec_int.size() -1 -i] = x[x.size() -1 -i];
		
		x = vec_int;
	}
	
	for (int i = 0; i < (x.size() - a.size()); i++)
	
	    { x[i] = vector<int>(1); x[i][0] = s; }
	
	
	//  Sort the vector
	
	vector<vector<int>> y (x.size());
	
	int N = x.size() / 2;
	
	for (int i = 1; i <= N; i <<= 1)
	{
		for (int j = 0; j < 2*N; j += 2*i)
		{
			for (int k = j, k1 = 0, k2 = 0; k < (j + 2*i); k++)
			{
				if ( (k1 < i) && (k2 < i) )
				{
					if (x[j+0 + k1][0] <= x[j+i + k2][0])
					
					      y[k] = x[j+0 + k1++];
					else  y[k] = x[j+i + k2++];
				}
				
				else
				{	if (k1 < i)
					
					      y[k] = x[j+0 + k1++];
					else  y[k] = x[j+i + k2++];
				}
			}
		}
		
		for (int k = 0; k < x.size(); k++) x[k] = y[k];
	}
	
	//  Unpad the vector
	
	for (int i = 0; i < a.size(); i++)
	
	    a[i] = x[(x.size() - a.size()) + i];
}



void Math::sort(int a[], int elements)
{
	//  Copy and pad the array
	
	int a_length = elements;
	
	int temp[a_length];
	
	int *x = temp;
	
	int x_length = a_length;
	
	for (int i = 0; i < a_length; i++)
	
	    x[i] = a[i];
	
	//  Find the smallest element
	
	int s = 0;
	
	for (int i = 0; i < a_length; i++)
	
	    if (a[i] < s) s = a[i];
	
	int *array1 = nullptr;
	
	if ((a_length & (a_length - 1)) != 0)
	{
		//  Pad the array to a power of two
		
		int length = a_length;
		
		int poweroftwo = 1;
		
		while (length != 0)
		{
			length >>= 1;
			
			poweroftwo <<= 1;
		}
		
		array1 = new int[poweroftwo];
		
		for (int i = 0; i < a_length; i++)
		
		    array1[poweroftwo -1 -i]
		       = x[a_length   -1 -i];
		
		x = array1;
		
		x_length = poweroftwo;
	}
	
	for (int i = 0; i < x_length - a_length; i++)
	
	    x[i] = s;
	
	
	//  Sort the array
	
	int y[x_length];
	
	int N = x_length / 2;
	
	for (int i = 1; i <= N; i <<= 1)
	{
		for (int j = 0; j < 2*N; j += 2*i)
		{
			for (int k = j, k1 = 0, k2 = 0; k < (j + 2*i); k++)
			{
				if ((k1 < i) && (k2 < i))
				{
					if (x[j + 0 + k1] <= x[j + i + k2])
					
					     y[k] = x[j + 0 + k1++];
					else y[k] = x[j + i + k2++];
				}
				
				else
				{	if (k1 < i)
					
					     y[k] = x[j + 0 + k1++];
					else y[k] = x[j + i + k2++];
				}
			}
		}
		
		for (int k = 0; k < x_length; k++)
		
		    x[k] = y[k];
	}
	
	//  Unpad the array
	
	for (int i = 0; i < a_length; i++)
	
	    a[i] = x[(x_length - a_length) + i];
	
	delete[] array1;
}


void Math::sort(int *a[], int elements)
{
	//  Copy and pad the array
	
	int **x;
	
	int a_length = elements;
	int x_length = elements;
	
	//  Find the smallest element
	
	int s = 0;
	
	for (int i = 0; i < elements; i++)
	{
		int int1 = a[i][0];
		
		if (int1 < s) s = int1;
	}
	
	int **array = nullptr;
	
	if ((x_length & (x_length -1)) != 0)
	{
		//  Pad the array to a power of two
		
		int length = x_length;
		
		int poweroftwo = 1;
		
		while (length != 0)
		
		    { length >>= 1; poweroftwo <<= 1; }
		
		array = new int*[poweroftwo];
		
		for (int i = 0; i < poweroftwo; i++) array[i] = nullptr;
		
		for (int i = 0; i < poweroftwo; i++)
		
		    array[i] = new int[2] { 0, 0 };
		
		for (int i = 0; i < x_length; i++)
		{
			array[poweroftwo -1 -i][0] = a[a_length -1 -i][0];
			array[poweroftwo -1 -i][1] = a[a_length -1 -i][1];
		}
		
		x = array;  x_length = poweroftwo;
	}
	
	else
	{	x = new int*[elements];
		
		for (int i = 0; i < elements; i++)
		
		    x[i] = new int[2] { a[i][0], a[i][1] };
	}
	
	
	//  Sort the array
	
	int **y = new int*[x_length];
	
	for (int i = 0; i < x_length; i++)
	
	    y[i] = new int[2] { 0, 0 };
	
	
	int N = x_length / 2;
	
	for (int i = 1; i <= N; i <<= 1)
	{
		for (int j = 0; j < 2*N; j += 2*i)
		{
			for (int k = j, k1 = 0, k2 = 0; k < (j + 2*i); k++)
			{
				if ( (k1 < i) && (k2 < i) )
				{
					if (x[j+0 + k1][0] <= x[j+i + k2][0])
					
					     { y[k][0] = x[j+0 + k1][0]; y[k][1] = x[j+0 + k1++][1]; }
					else { y[k][0] = x[j+i + k2][0]; y[k][1] = x[j+i + k2++][1]; }
				}
				
				else
				{	if (k1 < i)
					
					     { y[k][0] = x[j+0 + k1][0]; y[k][1] = x[j+0 + k1++][1]; }
					else { y[k][0] = x[j+i + k2][0]; y[k][1] = x[j+i + k2++][1]; }
				}
			}
		}
		
		for (int k = 0; k < x_length; k++) { x[k][0] = y[k][0];  x[k][1] = y[k][1]; }
	}
	
	
	//  Unpad the array
	
	for (int i = 0; i < a_length; i++)
	{
		a[i][0] = x[(x_length - a_length) + i][0];
		a[i][1] = x[(x_length - a_length) + i][1];
	}
	
	
	//  Delete the arrays
	
	for (int i = 0; i < x_length; i++)
	
	    { delete[] x[i]; x[i] = 0; }
	
	for (int i = 0; i < x_length; i++)
	
	    delete[] y[i];
	
	
	//  Delete the pointers to arrays
	
	//  (Don't write delete[] x, y)
	
	delete[] x;  delete[] y;
}


void Math::sort(Number a[], int elements)
{
	//  Copy and pad the array
	
	int a_length = elements;
	
	Number *x; Number temp[a_length]; x = temp;
	
	int x_length = a_length;
	
	for (int i = 0; i < a_length; i++)
	
	    x[i] = a[i];
	
	//  Find the smallest element
	
	Number s(0);
	
	for (int i = 0; i < a_length; i++)
	
	    if (a[i].less_than(s)) s = a[i];
	
	if ((a_length & (a_length - 1)) != 0)
	{
		//  Pad the array to a power of two
		
		int length = a_length;
		
		int poweroftwo = 1;
		
		while (length != 0)
		{
			length >>= 1;
			
			poweroftwo <<= 1;
		}
		
		Number array1[poweroftwo];
		
		for (int i = 0; i < a_length; i++)
		
		    array1[poweroftwo -1 -i]
		       = x[a_length   -1 -i];
		
		x = array1;
		
		x_length = poweroftwo;
	}
	
	for (int i = 0; i < x_length - a_length; i++)
	
	    x[i] = s;
	
	
	//  Sort the array
	
	Number y[x_length];
	
	int N = x_length / 2;
	
	for (int i = 1; i <= N; i <<= 1)
	{
		for (int j = 0; j < 2 * N; j += 2 * i)
		{
			for (int k = j, k1 = 0, k2 = 0; k < (j + 2 * i); k++)
			{
				if ((k1 < i) && (k2 < i))
				{
					if (x[j + 0 + k1].compare(x[j + i + k2]) <= 0)
					
					     y[k] = x[j + 0 + k1++];
					else y[k] = x[j + i + k2++];
				}
				
				else
				{	if (k1 < i)
					
					     y[k] = x[j + 0 + k1++];
					else y[k] = x[j + i + k2++];
				}
			}
		}
		
		for (int k = 0; k < x_length; k++)
		
		    x[k] = y[k];
	}
	
	//  Unpad the array
	
	for (int i = 0; i < a_length; i++)
	
	    a[i] = x[(x_length - a_length) + i];
}


vector<vector<int>> Math::sort_and_collate(vector<int> a)
{

	//  returns the sorted vector and the num-
	//  ber of times that each element occurs
	//  Duplicates of elements are removed.
	
	
	//  For example the array
	//
	//   { 3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8
	//     9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4 };
	//
	//  is sorted and collated to
	//
	//  [ 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4,
	//    5, 5, 5, 6, 6, 6, 7, 8, 8, 9, 9, 9 ]  or
	//
	//  { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 3 }, { 5, 3 },
	//  { 6, 3 }, { 7, 1 }, { 8, 2 }, { 9, 3 }
	
	
	Math::sort(a); // Sort the vector
	
	
	//  Count the number of unique elements
	
	int i, j, k;
	
	for (i = 1, j = 1; j < a.size(); j++)
	
	    if (a[j] != a[j - 1]) i++;
	
	int int_freq[i][2];
	
	int int_freq_length = i;
	
	//  Find the frequency of each element
	
	i = 0; j = 0; k = a[0];
	
	for (int s = 0; s < a.size(); s++)
	{
		int int1 = a[s];
		
		if (int1 == k) j++;
		
		else
		{	int_freq[i][0] = k;
			 int_freq[i][1] = j;
			
			i++; j = 1;
		}
		
		k = int1;
	}
	
	int_freq[int_freq_length - 1][0] = k;
	int_freq[int_freq_length - 1][1] = j;
	
	vector<vector<int>> int_freq1 (int_freq_length);
	
	for (int i = 0; i < int_freq1.size(); i++)
	
	    int_freq1[i] = vector<int>(2);
	
	for (int i = 0; i < int_freq_length; i++)
	{
		int_freq1[i][0] = int_freq[i][0];
		int_freq1[i][1] = int_freq[i][1];
	}
	
	return int_freq1;
}



vector<int*> Math::sort_and_collate(int a[], int elements)
{

	//  returns the sorted array and the num-
	//  ber of times that each element occurs
	//  Duplicates of elements are removed.
	
	
	//  For example the array
	//
	//   { 3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8
	//     9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4 };
	//
	//  is sorted and collated to
	//
	//  [ 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4,
	//    5, 5, 5, 6, 6, 6, 7, 8, 8, 9, 9, 9 ]  or
	//
	//  { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 3 }, { 5, 3 },
	//  { 6, 3 }, { 7, 1 }, { 8, 2 }, { 9, 3 }
	
	
	int a_length = elements;
	
	
	//  Make a copy of the array
	
	int temp_vector[a_length];
	
	for (int i = 0; i < a_length; i++)
	
	    temp_vector[i] = a[i];
	
	a = temp_vector;
	
	
	//  Sort the array
	
	sort(a, a_length);
	
	
	//  Count the number of unique elements
	
	int i, j, k;
	
	for (i = 1, j = 1; j < a_length; j++)
	
	    if (a[j] != a[j - 1]) i++;
	
	int int_freq[i][2];
	
	int int_freq_length = i;
	
	//  Find the frequency of each element
	
	i = 0; j = 0; k = a[0];
	
	for (int s = 0; s < a_length; s++)
	{
		int int1 = a[s];
		
		if (int1 == k) j++;
		
		else
		{	int_freq[i][0] = k;
			int_freq[i][1] = j;
			
			i++; j = 1;
		}
		
		k = int1;
	}
	
	int_freq[int_freq_length - 1][0] = k;
	int_freq[int_freq_length - 1][1] = j;
	
	vector<int*> int_freq1;
	
	for (int i = 0; i < int_freq_length; i++)
	
	    int_freq1.push_back(int_freq[i]);
	
	return int_freq1;
}


double Math::sqrt(double n)
{
	return root(n, 2);
}

double Math::square(double n)
{
	return n * n;
}


vector<int> Math::subtract(vector<int> minuend, vector<int> subtrahend)
{
	//  subtracts two int vectors
	
	vector<int> twos_comp_vec = twos_complement(subtrahend);
	
	vector<int> difference = add(minuend, twos_comp_vec);
	
	return difference;
}


int *Math::subtract(int minuend[], int subtrahend[], int elements)
{
	//  subtracts two int arrays
	
	int *twos_comp_array = twos_complement(subtrahend, elements);
	
	int *difference = add(minuend, twos_comp_array, elements);
	
	delete[] twos_comp_array;
	
	return difference;
}


double Math::tan(double x)
{
	return sin(x) / cos(x);
}

double Math::tanh(double x)
{
	return sinh(x) / cosh(x);
}


bool Math::test_bit(vector<int> a, long bit)
{
	int a_length = a.size();
	
	if (bit / 32 >= a_length) return false;
	
	return (a[a_length -1 - (int) (bit / 32)]
	
	    & (1 << (int) (bit % 32))) != 0 ? true : false;
}

bool Math::test_bit(int array[], int elements, long bit)
{
	int array_length = elements;
	
	if (bit / 32 >= array_length) return false;
	
	return (array[array_length -1 - (int) (bit / 32)]
	
	    & (1 << (int) (bit % 32))) != 0 ? true : false;
}


vector<int> Math::trim(vector<int> vec_int)
{

	//  removes the leading zeros from a vector
	//  If the elements are all zeros, it returns a one-int vector
	
	int vec_length = vec_int.size();
	
	if (vec_length == 1)
	{
		vector<int> vec_int1 (1);
		
		vec_int1[0] = vec_int[0];
		
		return vec_int1;
	}
	
	if (vec_length == 0)
	{
		vector<int> vec_int1 { 0 };
		
		return vec_int1;
	}
	
	//  Count the number of leading zeros
	
	int index = 0, zeros = 0;
	
	while ((index < vec_length - 1) && (vec_int[index] == 0))
	
	    { index++;  zeros++; }
	
	int size = vec_length - zeros;
	
	if (size == 0)  size = 1;
	
	vector<int> vec_int1 (size);
	
	int vec_length1 = vec_int1.size();
	
	for (int i = 0; i < vec_length1; i++)
	
	    vec_int1[vec_length1 -1 -i]
	  = vec_int [vec_length  -1 -i];
	
	return vec_int1;
}


int *Math::trim(int array[], int elements)
{
	//  removes the leading zeros from an array
	
	//  If the elements are all zeros, it returns a one-int array
	
	int array_length = elements;
	
	if (array_length == 1)
	{
		int *array1 = new int[1];
		
		array1[0] = array[0];
		
		return array1;
	}
	
	//  Count the number of leading zeros
	
	int index = 0, zeros = 0;
	
	while ((index < array_length - 1) && (array[index] == 0))
	
	    { index++;  zeros++; }
	
	int size = array_length - zeros;
	
	if (size == 0)  size = 1;
	
	int *array1 = new int[size];
	
	int array_length1 = size;
	
	for (int i = 0; i < array_length1; i++)
	
	    array1[array_length1 -1 -i]
	  = array [array_length  -1 -i];
	
	return array1;
}


vector<int> Math::trim(vector<int> vec_int, int zeros)
{
	//  removes the specified number of leading zeros
	
	int elements = vec_int.size();
	
	if (zeros >= elements)
	{
		string message = "zeros >= vector size";
		
		cout << message << endl; throw string(message);
	}
	
	int vec_length = elements;
	
	if (vec_length == 1)
	{
		vector<int> vec_int1 (1);
		
		vec_int1[0] = vec_int[0];
		
		return vec_int1;
	}
	
	int size = vec_length - zeros;
	
	vector<int> vec_int1 (size);
	
	int vec_length1 = size;
	
	for (int i = 0; i < vec_length1; i++)
	
	    vec_int1[vec_length1 -1 -i]
	  = vec_int [vec_length  -1 -i];
	
	return vec_int1;
}



int *Math::trim(int array[], int elements, int zeros)
{
	//  removes the specified number of leading zeros
	
	if (zeros >= elements)
	{
		string message = "zeros >= array size";
		
		cout << message << endl; throw string(message);
	}
	
	int array_length = elements;
	
	if (array_length == 1)
	{
		int *array1 = new int[1];
		
		array1[0] = array[0];
		
		return array1;
	}
	
	int size = array_length - zeros;
	
	int *array1 = new int[size];
	
	int array_length1 = size;
	
	for (int i = 0; i < array_length1; i++)
	
	    array1[array_length1 -1 -i]
	  = array [array_length  -1 -i];
	
	return array1;
}



vector<int> Math::twos_complement(vector<int> vec_int)
{
	//  multiplies a vector by -1
	
	//  Complement the bits and add 1
	
	int elements = vec_int.size();
	
	vector<int> vec_int1 (elements);
	
	for (int i = 0; i < elements; i++)
	
	    vec_int1[i] = ~vec_int[i];
	
	add_bit(vec_int1, 0);
	
	return vec_int1;
}


int *Math::twos_complement(int array[], int elements)
{
	//  multiplies an array by -1
	
	//  Complement the bits and add 1
	
	int *array1 = new int[elements];
	
	for (int i = 0; i < elements; i++)
	
	    array1[i] = ~array[i];
	
	add_bit(array1, elements, 0);
	
	return array1;
}
















const string Convert::base_16_separator = "0123456789abcdef";


const char Convert::int_to_base_64[64]
{
	'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P',
	'Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f',
	'g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v',
	'w','x','y','z','0','1','2','3','4','5','6','7','8','9','+','/',
};

const int Convert::base_64_to_int[128 - 5]
{
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 62, -1, -1, -1, 63,
	52, 53, 54, 55, 56, 57, 58, 59, 60, 61, -1, -1, -1, -1, -1, -1,
	-1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
	15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, -1, -1, -1, -1, -1,
	-1, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
	41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51
};



string Convert::char_array_to_string(char array[], int elements)
{
	string str;
	
	int array_length = elements;
	
	for (int i = 0; i < array_length; i++)
	
	    str + array[i];
	
	return str;
}


char  *Convert::string_to_char_array(string &str)
{
	int array_length = str.length();
	
	char *charray = new char[array_length];
	
	for (int i = 0; i < array_length; i++)
	
	    charray[i] = str[i];
	
	return charray;
}


string Convert::char_array_to_base_64(char array[], int elements)
{
	//  converts an array of octets to a string of sextets
	
	if ((elements % 3) != 0) throw string("array length is not divisible by 3"); 
	
	char *char64array = char_256_array_to_char_64_array(array, elements);
	
	int charray_length = elements;
	
	char charray[charray_length];
	
	static const char *int_to_base_64 = int_to_base_64;
	
	for (int i = 0; i < charray_length; i++)
	
	    charray[i] = int_to_base_64[char64array[i]];
	
	string str = string(charray);
	
	while ((str.length() % 4) != 0)  str += "=";
	
	string str1 (str);
	
	return str1;
}


char *Convert::base_64_to_char_array(string &str)
{
	//  converts a string of sextets to an array of octets
	
	if ((str.length() % 4) != 0) throw string("string length is not divisible by 4");
	
	if (!Number::is_base_64(str))  throw string("string is not in base 64");
	
	int array_length = str.length();
	
	char *array = new char[array_length];
	
	static const int *base64toint = base_64_to_int;
	
	for (int i = 0; i < array_length; i++)
	
	    array[i] = base64toint[str[i]];
	
	array = char_64_array_to_char_256_array(array, str.length());
	
	return array;
}


char *Convert::char_64_array_to_char_256_array(char array[], int elements)
{
	//  converts a 6-bit-char array to an 8-bit-char array
	
	if ((elements % 4) != 0) throw string("array length not divisible by 4"); 
	
	int array_length = elements;
	
	for (int i = 0; i < array_length -2; i++)
	{
		if ((array[i] & 0xc0) != 0)
		{
			string message = "array is not in base 64";
			
			throw string(message);
		}
	}
	
	//  Two equal chars == 1 char  == 1 x 6 + 2 bits modulo 3
	//  One equal char  == 2 chars == 2 x 6 + 4 bits modulo 3
	//   No equal char  == 3 chars == 3 x 6 + 6 bits modulo 3
	
	//  Count the equals chars
	
	int equals = 0;
	
	for (int i = 0; i < 2; i++)
	
	    if (      array_length -1 -i   >= 0)
	    if (array[array_length -1 -i] == -1) equals++;
	
	char char64array[array_length];
	
	for (int i = 0; i < array_length; i++)
	
	    char64array[i] = array[i];
	
	int char256_length = array_length * 3/4 - equals;
	
	char *char256array = new char[char256_length];
	
	//  Example  Convert 4 sextets into 3 octets
	//
	//  XX101101 XX010101 XX111111 XX001100
	//    101101   01|0101  1111|11  001100
	//    10110101    01011111   11001100
	
	for (int i = 0, j = 0; i < (4 + char256_length/3); i++, j++)
	{
		char carray[3]; for (int k = 0; k < 3; k++) carray[k] = 0;
		
		if ((4*j + 0) < array_length)
		
		    carray[0] = (char) (char64array[4*j + 0] << 2);
		
		if ((4*j + 1) < array_length)
		{
			carray[0] +=      ((char64array[4*j + 1] >> 4) & 0x3f);
			carray[1] = (char) (char64array[4*j + 1] << 4);
		}
		
		if ((4*j + 2) < array_length)
		{
			carray[1] +=      ((char64array[4*j + 2] >> 2) & 0x0f);
			carray[2] = (char) (char64array[4*j + 2] << 6);
		}
		
		if ((4*j + 3) < array_length)
		
		    carray[2] += ((char64array[4*j + 3] >> 0) & 0x3f);
		
		for (int k = 0; k < 3; k++)
		
		    if ((3*i + k) < char256_length)
		
			char256array[3*i + k] += carray[k];
	}
	
	return char256array;
}


char *Convert::char_256_array_to_char_64_array(char array[], int elements)
{
	//  converts an 8-bit-char array to a 6-bit-char array
	
	if ((elements % 3) != 0) throw string("array is not divisible by 3");
	
	int array_length = elements;
	
	char char256array[array_length];
	
	for (int i = 0; i < array_length; i++)
	
	    char256array[i] = array[i];
	
	int r = array_length % 3;
	
	int char64_length = (array_length - r) * 4/3;
	
	if      (r == 1) char64_length += 2; // 8 bits == 1x6 + 2 bits
	else if (r == 2) char64_length += 3; // 16 bits == 2 x 6 + 4 bits
	
	char *char64array = new char[char64_length];
	
	
	//  Example  Convert 3 octets into 4 sextets (3*8 == 4*6)
	//
	//    10110101    01011111     11001100
	//    101101   01|0101  1111|11  001100
	//  XX101101 XX010101 XX111111 XX001100
	
	//  char256  0 1 2    3 4 5    6 7 8       9 10 11 ...
	//  char64   0 1 2 3  4 5 6 7  8 9 10 11  12 13 14 15 ...
	
	for (int i = 0, j = 0; i < 1 + char64_length / 4; i++, j++)
	{
		char carray[4]; for (int k = 0; k < 4; k++) carray[k] = 0;
		
		if ((3*j + 0) < array_length)
		{
			carray[0] =  (char)((char256array[3*j + 0] >> 2) & 0x3f);
			carray[1] =  (char)((char256array[3*j + 0] << 4) & 0x30);
		}
		
		if ((3*j + 1) < array_length)
		{
			carray[1] +=       ((char256array[3*j + 1] >> 4) & 0x0f);
			carray[2] =  (char)((char256array[3*j + 1] << 2) & 0x3c);
		}
		
		if ((3*j + 2) < array_length)
		{
			carray[2] += (char)((char256array[3*j + 2] >> 6) & 0x03);
			carray[3] =  (char)((char256array[3*j + 2] << 0) & 0x3f);
		}
		
		for (int k = 0; k < 4; k++)
		
		    if ((4*i + k) < char64_length)
		
			char64array[4*i + k] += carray[k];
	}
	
	return char64array;
}


char *Convert::char_16_array_to_char_256_array(char array[], int elements)
{
	//  converts a 4-bit-char array to an
	//
	//  8-bit-char array and halves the size of the array
	
	int array_length = elements;
	
	int char16_length = ((array_length & 1) == 0) ?
	
	    array_length : 1 + array_length;
	
	char char16array[char16_length];
	
	for (int i = 0; i < char16_length; i++) char16array[i] = 0;
	
	for (int i = 0; i < array_length; i++)
	
	    char16array[i + char16_length - array_length] = array[i];
	
	int char256_length = char16_length / 2;
	
	char *char256array = new char[char256_length];
	
	for (int i = 0; i < char256_length; i++) char256array[i] = 0;
	
	for (int i = 0; i < char16_length; i++)
	{
		if (char16array[i] >= 'a')
		{
			char16array[i] -= 'a';
			char16array[i] += 10;
		}
		
		else char16array[i]-='0';
	}
	
	for (int i = 0; i < char256_length; i++)
	{
		int temp = 0;
		temp += char16array[2*i];
		temp &= 0xf;
		
		char256array[i] += temp;
		char256array[i] <<= 4;
		
		temp = 0;
		temp += char16array[2*i + 1];
		temp &= 0xf;
		
		char256array[i] += temp;
	}
	
	return char256array;
}



vector<char> Convert::char_256_vector_to_char_16_vector(vector<char> vec_char)
{
	//  converts an 8-bit-char vector to a
	//
	//  4-bit-char vector and doubles the size of the vector
	
	int char_256_length = vec_char.size();
	
	vector<char> char_256_vec (char_256_length);
	
	for (int i = 0; i < char_256_length; i++)
	
	    char_256_vec[i] = vec_char[i];
	
	int char_16_length = char_256_length * 2;
	
	vector<char> char_16_vec (char_16_length);
	
	for (int i = 0; i < char_16_length; i++) char_16_vec[i] = 0;
	
	for (int i = 0; i < char_256_length; i++)
	{
		char_16_vec[2*i + 1] += char_256_vec[i];
		char_16_vec[2*i + 1] &= 0xf;
		
		char_256_vec[i] >>= 4;
		char_256_vec[i] &= 0xf;
		
		char_16_vec[2*i + 0] += char_256_vec[i];
		char_16_vec[2*i + 0] &= 0xf;
	}
	
	//  character values of the integer values
	
	char hex_char_array[16] =
	
	    { '0', '1', '2', '3', '4', '5', '6', '7',
	      '8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };
	
	//  Convert from int values 0 to 15 to char values '0' to 'f'
	
	for (int i = 0; i < char_16_length; i++)
	
	    char_16_vec[i] = hex_char_array[char_16_vec[i]];
	
	return char_16_vec;
}



char *Convert::char_256_array_to_char_16_array(char array[], int elements)
{
	//  converts an 8-bit-char array to a
	//
	//  4-bit-char array and doubles the size of the array
	
	int array_length = elements;
	
	char char256array[array_length];
	
	for (int i = 0; i < array_length; i++)
	
	    char256array[i] = array[i];
	
	int char16_length = array_length * 2;
	
	char *char16array = new char[char16_length];
	
	for (int i = 0; i < char16_length; i++) char16array[i] = 0;
	
	for (int i = 0; i < array_length; i++)
	{
		char16array[2*i + 1] += char256array[i];
		char16array[2*i + 1] &= 0xf;
		
		char256array[i] >>= 4;
		char256array[i] &= 0xf;
		
		char16array[2*i + 0] += char256array[i];
		char16array[2*i + 0] &= 0xf;
	}
	
	//  character values of the integer values
	
	char hex_char_array[16] =
	
	    { '0', '1', '2', '3', '4', '5', '6', '7',
	      '8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };
	
	//  Convert from int values 0 to 15 to char values '0' to 'f'
	
	for (int i = 0; i < char16_length; i++)
	
	    char16array[i] = hex_char_array[char16array[i]];
	
	return char16array;
}



int *Convert::char_array_to_int_array(char array[], int elements)
{
	//  converts a char array to an int array
	
	int array_length = (elements + 3) / 4;
	
	int *intarray = new int[array_length];
	
	for (int i = 0; i < array_length; i++) intarray[i] = 0;
	
	for (int i = 0; i < array_length; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			int temp = 0;
			
			temp += array[4*i + j];
			
			temp &= 0xff;
			
			intarray[i] += temp;
			
			if (j < 3) intarray[i] <<= 8;
		}
	}
	
	return intarray;
}



vector<char> Convert::int_vector_to_char_vector(vector<int> vec_int)
{
	//  converts an int vector to a char vector
	
	int vec_int_length = vec_int.size();
	
	vector<int> vec_int1 (vec_int_length);
	
	for (int i = 0; i < vec_int_length; i++)
	
	    vec_int1[i] = vec_int[i];
	
	int vec_char_length = vec_int_length * 4;
	
	vector<char> vec_char (vec_char_length);
	
	for (int i = 0; i < vec_char_length; i++) vec_char[i] = 0;
	
	for (int i = 0; i < vec_int_length; i++)
	{
		vec_char[4*i + 3] += vec_int[i];  vec_int[i] >>= 8;  vec_int[i] &= 0x00ffffff;
		vec_char[4*i + 2] += vec_int[i];  vec_int[i] >>= 8;
		vec_char[4*i + 1] += vec_int[i];  vec_int[i] >>= 8;
		vec_char[4*i + 0] += vec_int[i];
	}
	
	return vec_char;
}



char *Convert::int_array_to_char_array(int array[], int elements)
{
	//  converts an int array to a char array
	
	int array_length = elements;
	
	int intarray[array_length];
	
	for (int i = 0; i < array_length; i++)
	
	    intarray[i] = array[i];
	
	int charray_length = array_length * 4;
	
	char *charray = new char[charray_length];
	
	for (int i = 0; i < charray_length; i++) charray[i] = 0;
	
	for (int i = 0; i < array_length; i++)
	{
		charray[4*i + 3] += intarray[i];  intarray[i] >>= 8;  intarray[i] &= 0x00ffffff;
		charray[4*i + 2] += intarray[i];  intarray[i] >>= 8;
		charray[4*i + 1] += intarray[i];  intarray[i] >>= 8;
		charray[4*i + 0] += intarray[i];
	}
	
	return charray;
}


int *Convert::double_array_to_int_array(double array[], int elements)
{
	int array_length = elements;
	
	int *intarray = new int[array_length];
	
	for (int i = 0; i < array_length; i++)
	
	    intarray[i] = (int) array[i];
	
	return intarray;
}


double *Convert::double_array_to_complex_double_array(double array[], int elements)
{
	int array_length = elements;
	
	int array1_length = array_length * 2;
	
	double *array1 = new double[array1_length];
	
	for (int i = 0; i < array1_length; i++) array1[i] = 0;
	
	for (int i = 0; i < array_length; i++)
	{
		array1[2*i + 0] = array[i];
		array1[2*i + 1] = 0.0D;
	}
	
	return array1;
}


double *Convert::complex_double_array_to_double_array(double array[], int elements)
{
	int array_length = elements;
	
	int array1_length = array_length / 2;
	
	double *array1 = new double[array1_length];
	
	for (int i = 0; i < array1_length; i++) array1[i] = 0;
	
	for (int i = 0; i < array_length / 2; i++)
	
	    array1[i] = array[2*i + 0];
	
	return array1;
}


double *Convert::double_array(double array1[], int elements)
{

	//  doubles the length of a double array
	
	double *d_array;
	
	int array_length1 = elements;
	
	int d_length = 2 * array_length1;
	
	d_array = new double[d_length];
	
	for (int i = 0; i < d_length; i++) d_array[i] = 0;
	
	for (int i = 0; i < array_length1; i++)
	
	    d_array[i] = array1[i];
	
	return d_array;
}



//  The compiler says that ISO C++ forbids converting from array to vector or vice versa
//  (to prevent the user from accidentally assigning one type of object to another type).
//  This means that the user has to do the conversion explicitly because the compiler
//  won't allow it to be done implicitly.


vector<int> Convert::int_array_to_int_vector(int array[], int elements)
{
	vector<int> vec_int (elements);
	
	for (int i = 0; i < elements; i++)
	
	    vec_int[i] = array[i];
	
	return vec_int;
}


int* Convert::int_vector_to_int_array(vector<int> vec_int)
{
	int *array = new int[vec_int.size()];
	
	for (int i = 0; i < vec_int.size(); i++)
	
	    array[i] = vec_int[i];
	
	return array;
}






int main()
{





}





