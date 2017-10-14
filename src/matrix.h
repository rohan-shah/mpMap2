#ifndef MATRIX_HEADER_GUARD
#define MATRIX_HEADER_GUARD
struct integerMatrix
{
public:
	integerMatrix(const integerMatrix& mat)
		: data(mat.data), nRow(mat.nRow), nCol(mat.nCol)
	{}
	integerMatrix()
		: data(NULL), nRow(0), nCol(0)
	{}
	integerMatrix(int* data, std::size_t nRow, std::size_t nCol)
		: data(data), nRow(nRow), nCol(nCol)
	{}
	integerMatrix(Rcpp::IntegerMatrix mat)
		: data(&(mat[0])), nRow(mat.nrow()), nCol(mat.ncol())
	{}
	integerMatrix& operator=(Rcpp::IntegerMatrix& mat)
	{
		data = &(mat[0]);
		nRow = mat.nrow();
		nCol = mat.ncol();
		return *this;
	}
	int& operator()(std::size_t row, std::size_t column)
	{
		return data[column*nRow + row];
	}
	int* data;
	std::size_t nRow, nCol;
};
struct numericMatrix
{
public:
	numericMatrix(const numericMatrix& mat)
		: data(mat.data), nRow(mat.nRow), nCol(mat.nCol)
	{}
	numericMatrix()
		: data(NULL), nRow(-1), nCol(-1)
	{}
	numericMatrix(double* data, std::size_t nRow, std::size_t nCol)
		: data(data), nRow(nRow), nCol(nCol)
	{}
	numericMatrix& operator=(Rcpp::NumericMatrix& mat)
	{
		data = &(mat[0]);
		nRow = mat.nrow();
		nCol = mat.ncol();
		return *this;
	}
	numericMatrix(Rcpp::NumericMatrix mat)
		: data(&(mat[0])), nRow(mat.nrow()), nCol(mat.ncol())
	{}
	double& operator()(std::size_t row, std::size_t column)
	{
		return data[column*nRow + row];
	}
	double* data;
	std::size_t nRow, nCol;
};
#endif
