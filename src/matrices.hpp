#ifndef MATRICES_HEADER_GUARD
#define MATRICES_HEADER_GUARD
#include <vector>
template<typename T> class rowMajorMatrix
{
public:
	rowMajorMatrix()
		:nRows(0), nColumns(0)
	{}
	rowMajorMatrix(int nRows, int nColumns)
		: nRows(nRows), nColumns(nColumns), data(nRows*nColumns)
	{}
	rowMajorMatrix(int nRows, int nColumns, T value)
		: nRows(nRows), nColumns(nColumns), data(nRows*nColumns, value)
	{}
	rowMajorMatrix(rowMajorMatrix<T>&& other)
		:nRows(other.nRows), nColumns(other.nColumns), data(std::move(other.data))
	{}
	rowMajorMatrix<T>& operator=(rowMajorMatrix<T>&& other)
	{
		nRows = other.nRows;
		nColumns = other.nColumns;
		data.swap(other.data);
		return *this;
	}
	rowMajorMatrix<T> copy() const
	{
		rowMajorMatrix<T> retVal;
		retVal.nRows = nRows;
		retVal.nColumns = nColumns;
		retVal.data = data;
		return retVal;
	}
	typename std::vector<T>::iterator iterator(int row, int column)
	{
		return data.begin() + column + row*nColumns;
	}
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
	rowMajorMatrix(const rowMajorMatrix<T>& other);
	rowMajorMatrix<T>& operator=(const rowMajorMatrix<T>& other);
	int nRows, nColumns;
	std::vector<T> data;
};
template<typename T> class xMajorMatrix
{
public:
	xMajorMatrix()
		:sizeX(0), sizeY(0), sizeZ(0)
	{}
	xMajorMatrix(int sizeX, int sizeY, int sizeZ)
		: sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), data(sizeX*sizeY*sizeZ)
	{}
	xMajorMatrix(int sizeX, int sizeY, int sizeZ, T value)
		: sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), data(sizeX*sizeY*sizeZ, value)
	{}
	xMajorMatrix(xMajorMatrix<T>&& other)
		:sizeX(other.sizeX), sizeY(other.sizeY), sizeZ(other.sizeZ), data(std::move(other.data))
	{}
	xMajorMatrix<T>& operator=(xMajorMatrix&& other)
	{
		sizeX = other.sizeX;
		sizeY = other.sizeY;
		sizeZ = other.sizeZ;
		data.swap(other.data);
	}
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
	xMajorMatrix(const xMajorMatrix<T>& other);
	xMajorMatrix<T>& operator=(const xMajorMatrix<T>& other);
	int sizeX, sizeY, sizeZ;
	std::vector<T> data;
};
#endif
