#ifndef MATRICES_HEADER_GUARD
#define MATRICES_HEADER_GUARD
template<int n> struct array2
{
public:
	array2()
	{}
	double values[n][n];
};
template<typename T> class rowMajorMatrix
{
public:
	rowMajorMatrix(int nRows, int nColumns)
		: nRows(nRows), nColumns(nColumns), data(nRows*nColumns)
	{}
	rowMajorMatrix(int nRows, int nColumns, T value)
		: nRows(nRows), nColumns(nColumns), data(nRows*nColumns, value)
	{}
	typename std::vector<T>::reference operator()(int row, int column)
	{
		return data[column + row*nColumns];
	}
	typename std::vector<T>::const_reference operator()(int row, int column) const
	{
		return data[column + row*nColumns];
	}
	int getNRows() const
	{
		return nRows;
	}
	int getNColumns() const
	{
		return nColumns;
	}
	void swap(rowMajorMatrix<T>& other)
	{
		nRows = other.nRows;
		nColumns = other.nColumns;
		data.swap(other.data);
	}
private:
	int nRows, nColumns;
	std::vector<T> data;
};
template<typename T> class xMajorMatrix
{
public:
	xMajorMatrix(int sizeX, int sizeY, int sizeZ)
		: sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), data(sizeX*sizeY*sizeZ)
	{}
	xMajorMatrix(int sizeX, int sizeY, int sizeZ, T value)
		: sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), data(sizeX*sizeY*sizeZ, value)
	{}
	T& operator()(int i, int j, int k)
	{
		return data[i + j*sizeX + k * sizeX * sizeY];
	}
	int getSizeX()
	{
		return sizeX;
	}
	int getSizeY()
	{
		return sizeY;
	}
	int getSizeZ()
	{
		return sizeZ;
	}
	void swap(xMajorMatrix<T>& other)
	{
		sizeX = other.sizeX;
		sizeY = other.sizeY;
		sizeZ = other.sizeZ;
		data.swap(other.data);
	}
private:
	int sizeX, sizeY, sizeZ;
	std::vector<T> data;
};
#endif