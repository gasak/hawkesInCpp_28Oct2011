//
//  Matrix.h
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/11/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#include <iostream>
#include<new>
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include <iostream>


#ifndef hawkesInCpp_Matrix_h
#define hawkesInCpp_Matrix_h



using namespace std;

template<class T>
class Matrix
{
public:
    
	//1.Constructors
    
    
	Matrix(int Row_Size,int Column_Size);
    
	Matrix()
	{
		m_RowSize=0;
		m_ColumnSize=0;
	}
    
    
	Matrix(const Matrix<T> &);
    
	//1.end of Constructor declaration
    
	//2.destructor
    
	~Matrix(void)
	{
        
		for(int i=0;i<m_RowSize;i++)
			delete[] m_Data[i];
		delete [] m_Data;
	}
    
	//.end of destructor declaration
    
    
    
	//.4.Input Elements Function
    
	void Set_Element(int Row_,int Column,const T& Value);
    
	//end Input element Function declaration
    
	//5.Output Elements Function
    
	T Get_Element(int Row ,int Column) const ;
    
	//end Output Elements Function declaration
    
    
    
	//6.OPERATOR OVERLOADING
    
	//6.17 Unary Operator !
    
	bool operator!() const;
    
	//end unary operator declaration
    
	//6.16 Assigment Operator
    
	const Matrix<T> & operator=(const Matrix<T> &);
    
	const Matrix<T> & operator=(const T&);
    
	//end Assignment Operator
    
	//6.14 and 6.15 equality /non equality operators
    
	bool operator==(const Matrix<T> &) const;
    
	bool operator==(const T & ) const;
    
	bool operator!=(const Matrix<T> & second) const
    {
		return (!(*this==second));
    }
    
	bool operator!=(const T & number)
    {
		return !(*this==number);
    }
    
    
	//end equality/non equality operators declaration
    
	//6.13 data array operators
    
	T &operator() (int Row,int Column);
    
	//end data array operators declaration
    
	//6.3 and 6.1 + operator
    
	Matrix<T> operator+(const T &);
    
	Matrix<T> operator+(const Matrix<T> &);
    
	//end + operator declaration
    
	//6.4 += operator
    
	const Matrix<T> &operator+=(const T&);
    
	const Matrix<T> & operator+=(const Matrix<T> &);
    
	//end += operator declaration
    
	//6.5 and 6.7 - operator 
    
	Matrix<T> operator-(const T&);
    
	Matrix<T>  operator - (const Matrix<T> &);
    
	//end - operator declaration
    
	//6.8 -= operator 
	const Matrix<T> & operator -=(const T & );
    
	const Matrix<T> & operator -=(Matrix<T> & );
    
	//end -= operator declaration
    
	//6.9 Multiply Operator *
    
	Matrix<T>  operator*(const T&);
    
	Matrix<T>  operator* (const Matrix<T> &);
    
	//end Multiply operator Declaration
    
	//6.18 *= operator
    
	const Matrix<T> & operator*=(const T&);
    
	const Matrix<T> & operator*=(const Matrix<T> &);
    
	//END OPERATOR OVERLOADING
    
    
	void Transposed();
    
	Matrix<T>  GetTransposed() const;
    
	bool Inverse();
    
	Matrix<T> GetInversed() const ;
    
	double Determinant(double **a,int n) const;
	double Det() const;
    
	void Set_Size(int RowSize,int ColumnSize);
    
	const int& Get_ColumnSize() const;
    
	const int & Get_RowSize() const;
    
	void print() const;
    
	Matrix<T> getKthRow(int k);
    
	void setKthRow(int k, Matrix<T> & rowMatrix);
	double getRowNorm2(int k);
    
public:
    
	//3.Data declaration
    
	int m_RowSize;
	int m_ColumnSize;
	T**m_Data;
    
	//end Data declaration
    
};


template<class T>
double Matrix<T>::Determinant(double **a,int n) const{
	int i,j,j1,j2 ;                    // general loop and matrix subscripts
	double det = 0 ;                   // init determinant
	double **m = NULL ;                // pointer to pointers to implement 2d
	// square array
    
	if (n < 1)    {   }                // error condition, should never get here
    
	else if (n == 1) {                 // should not get here
		det = a[0][0] ;
	}
    
	else if (n == 2)  {                // basic 2X2 sub-matrix determinate
		// definition. When n==2, this ends the
		det = a[0][0] * a[1][1] - a[1][0] * a[0][1] ;// the recursion series
	}
    
    
	// recursion continues, solve next sub-matrix
	else {                             // solve the next minor by building a
		// sub matrix
		det = 0 ;                      // initialize determinant of sub-matrix
        
		// for each column in sub-matrix
		for (j1 = 0 ; j1 < n ; j1++) {
			// get space for the pointer list
			m = (double **) malloc((n-1)* sizeof(double *)) ;
            
			for (i = 0 ; i < n-1 ; i++)
				m[i] = (double *) malloc((n-1)* sizeof(double)) ;
			//     i[0][1][2][3]  first malloc
			//  m -> +  +  +  +   space for 4 pointers
			//       |  |  |  |          j  second malloc
			//       |  |  |  +-> _ _ _ [0] pointers to
			//       |  |  +----> _ _ _ [1] and memory for
			//       |  +-------> _ a _ [2] 4 doubles
			//       +----------> _ _ _ [3]
			//
			//                   a[1][2]
			// build sub-matrix with minor elements excluded
			for (i = 1 ; i < n ; i++) {
				j2 = 0 ;               // start at first sum-matrix column position
				// loop to copy source matrix less one column
				for (j = 0 ; j < n ; j++) {
					if (j == j1) continue ; // don't copy the minor column element
                    
					m[i-1][j2] = a[i][j] ;  // copy source element into new sub-matrix
					// i-1 because new sub-matrix is one row
					// (and column) smaller with excluded minors
					j2++ ;                  // move to next sub-matrix column position
				}
			}
            
			det += pow(-1.0,1.0 + j1 + 1.0) * a[0][j1] * Determinant(m,n-1) ;
			// sum x raised to y power
			// recursively get determinant of next
			// sub-matrix which is now one
			// row & column smaller
            
			for (i = 0 ; i < n-1 ; i++) free(m[i]) ;// free the storage allocated to
			// to this minor's set of pointers
			free(m) ;                       // free the storage for the original
			// pointer to pointer
		}
	}
	return(det) ;
}




template<class T>
double Matrix<T>::Det() const
{
	if(m_RowSize!=m_ColumnSize)
		exit(0);
	int n=m_RowSize;
	double **a=new double *[n];
	for(int i=0;i<n;i++)
		a[i]=new double [n];
    
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			a[i][j]=(double)m_Data[i][j];
	double det= Determinant(a,n);
    
	//Determinant(a,n);
	return det;
    
    
}


template<class T>
Matrix<T>::Matrix(int Row_Size,int Column_Size)
:m_RowSize(Row_Size),m_ColumnSize(Column_Size)
{
	m_Data=new T *[m_RowSize];
	for(int i=0;i<m_RowSize;i++)
	{
		m_Data[i]=new T [m_ColumnSize];
	}
    
}
//Set Size of the array with this Constructor
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
Matrix<T>::Matrix(const Matrix<T> & second)
{
	m_RowSize=second.Get_RowSize();
	m_ColumnSize=second.Get_ColumnSize();
	m_Data=new T *[m_RowSize];
	for(int i=0;i<m_RowSize;i++)
	{
		m_Data[i]=new T [m_ColumnSize];
	}
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<m_ColumnSize;j++)
			m_Data[i][j]=second.Get_Element(i,j);
}
//Copy a Matrix to the new Initialized one by using this constructor
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
void Matrix<T>::Set_Element(int Row,int Column,const T& Value)
{
	if(m_RowSize>Row && m_ColumnSize>Column)
		m_Data[Row][Column]=Value;
	//else exit(0);
}
//Set Element i,j of the matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
T Matrix<T>::Get_Element(int Row, int Column) const
{if(m_RowSize>Row && m_ColumnSize>Column)
	return m_Data[Row][Column];
else exit(0);
}
//Get Element i,j of the Matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
bool Matrix<T>::operator !() const
{
	return (m_RowSize==0);
}
//See weather the matrix is initialized or not
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
const Matrix<T> & Matrix<T>::operator =(const Matrix<T> & second)
{
	for(int i=0;i<m_RowSize;i++)
		delete[] m_Data[i];
	delete[] m_Data;
	m_RowSize=second.Get_RowSize();
	m_ColumnSize=second.Get_ColumnSize();
	m_Data=new T* [m_RowSize];
    
	for(int i=0;i<m_RowSize;i++)
		m_Data[i]=new T[m_ColumnSize];
    
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<m_ColumnSize;j++)
			m_Data[i][j]=second.Get_Element(i,j);
    
	return *this;
}
//Copy one Matrix to the other
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
bool Matrix<T>::operator ==(const Matrix<T> & second) const
{
	if(m_RowSize!=second.Get_RowSize() || m_ColumnSize!=second.Get_ColumnSize())
		return false;
	else
	{
		for(int i=0;i<m_RowSize;i++)
		{
			for(int j=0;j<m_ColumnSize;j++)
			{
				if(m_Data[i][j]!=second.Get_Element(i,j))
					return false;
			}
		}
	}
	return true;
}
//See Weather two Matrixes are equal
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
bool Matrix<T>::operator==(const T & number) const
{
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<m_ColumnSize;j++)
			if(m_Data[i][j]!=number)
				return false;
	return true;
}
//See Weather all elements of the matrix are equal to the specified number
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
T & Matrix<T>::operator ()(int Row,int Column)
{
	if(m_RowSize>Row && m_ColumnSize>Column)
		return m_Data[Row][Column];
	else exit(0);
}
//Get access to the element i,j
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
const Matrix<T> & Matrix<T>::operator +=(const T& number)
{
	for(int i=0;i<m_RowSize;i++)
	{
		for(int j=0;j<m_ColumnSize;j++)
		{
			m_Data[i][j]+=number;
		}
	}
    
	return *this;
}
//add the specified number to all the matrix elements
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
const Matrix<T> & Matrix<T>::operator+=(const Matrix<T> & second)
{
	if(m_RowSize!=second.Get_RowSize() || m_ColumnSize!=second.Get_ColumnSize())
        
		exit(0);
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<m_ColumnSize;j++)
		{
			m_Data[i][j]+=second.Get_Element(i,j);
		}
	return *this;
}
//Add another matrix to the Matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
Matrix<T>  Matrix<T>::operator +(const Matrix<T> & second)
{
	Matrix<T> Temp(*this);
	if(m_RowSize!=second.Get_RowSize() || m_ColumnSize!=second.Get_ColumnSize())
        
		exit(0);
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<m_ColumnSize;j++)
			Temp+=second;
	return Temp;
}
//Add 2 matrixes
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
Matrix<T> Matrix<T>::operator+(const T& number)
{
	Matrix<T> Temp(m_RowSize,m_ColumnSize);
    
	for(int i=0;i<m_RowSize;i++)
        
		for(int j=0;j<m_ColumnSize;j++)
            
			Temp(i,j)=Get_Element(i,j)+number;
    
	return Temp;
}
//Add a number to the matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
Matrix<T> Matrix<T>::operator -(const Matrix<T> & second)
{
	Matrix<T> Temp(second);
	if(m_RowSize!=second.Get_RowSize() || m_ColumnSize!=second.Get_ColumnSize())
		exit(0);
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<m_ColumnSize;j++)
			Temp(i,j)=-Temp(i,j);
	Temp+=*this;
	return Temp;
}
//Subtract a number from  matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
const Matrix<T>& Matrix<T>::operator -=(const T & number)
{
    
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<m_ColumnSize;j++)
			m_Data[i][j]-=number;
    
	return *this;
}
//Subtract a number from the Matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
const Matrix<T>  &Matrix<T>::operator -=(Matrix<T> & second)
{
	*this=(*this)-second;
	return *this;
}
//Subtract another Matrix from the Matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
Matrix<T>  Matrix<T>::operator *(const T & number)
{
	Matrix<T> Temp(*this);
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<m_ColumnSize;j++)
			Temp.Set_Element(i,j,Temp.Get_Element(i,j)*number);
	return Temp;
}
//Multiply a number to the matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
Matrix<T>  Matrix<T>::operator *(const Matrix<T> & second)
{
	if(m_ColumnSize!=second.Get_RowSize() || m_RowSize!=second.Get_ColumnSize())
		exit(0);
	Matrix<T> Temp(m_RowSize,second.Get_ColumnSize());
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<second.Get_ColumnSize();j++)
		{
			Temp(i,j)=0;
			for(int k=0;k<m_ColumnSize;k++)
				Temp(i,j)+=(m_Data[i][k]*second.Get_Element(k,j));
		}
	return Temp;
}
//Multiply 2 Matrixes
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>

const Matrix<T> &Matrix<T>::operator *=(const T &number)
{
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<m_ColumnSize;j++)
			m_Data[i][j]*=number;
	return *this;
}
//Multiply a number to the matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
const Matrix<T> & Matrix<T>::operator *=(const Matrix<T> & second)
{
	return ((*this)*second);
}
//Multiply a matrix to the matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
void Matrix<T>::Transposed()
{
	Matrix<T> Temp(m_ColumnSize,m_RowSize);
	for(int i=0;i<m_ColumnSize;i++)
		for(int j=0;j<m_RowSize;j++)
			Temp(i,j)=m_Data[j][i];
	*this=Temp;
}
//Make the matrix transpose
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>

Matrix<T>  Matrix<T>::GetTransposed() const
{
    
    
	Matrix<T> Temp(m_ColumnSize,m_RowSize);
	for(int i=0;i<m_ColumnSize;i++)
		for(int j=0;j<m_RowSize;j++)
			Temp.Set_Element(i,j,m_Data[j][i]);
    
	return Temp;
}
//Get transosed of the matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>

bool Matrix<T>::Inverse()
{   //This function uses gausian method to invert the function
	//This code might not be so clear because this function is the only one I have written it in code spagety !
	if(m_RowSize!=m_ColumnSize)
		exit(0);
	if(Det()==0)
		return false;
    
	Matrix<double> c(m_RowSize,m_ColumnSize*2);
    
    
    
	Matrix<double> I(m_RowSize,m_ColumnSize);
	int i,j,k,p=0;
    
	Matrix<double> tempar(2*m_RowSize,1);
	int cntRow=0;
	int r1=0,r2=0;
	Matrix<double> rowar(m_RowSize,2);
	double cii;
	int size=m_RowSize;
	for( i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			c(i,j)=(double)m_Data[i][j];
            
            
		}
	}
	for(i=0;i<size;i++)
		for(j=0;j<size;j++)
		{
			if(i!=j)
                
			{I(i,j)=0;
			}
			else if (i==j)
			{I(i,j)=1;
			}
		}
	for(i=0;i<size;i++)
		for(j=size;j<2*size;j++)
		{
			c(i,j)=I(i,j-size);
		}
	cntRow=0;
	for (i=0;i<size;i++)
	{
		r1=i;
		cii=c(i,i);
		for(j=i+1;j<size;j++)
		{
			if(c(j,i)>cii)
			{
				cii=c(j,i);
				r2=j;
			}
		}
		if(r2!=r1)
		{
			rowar(cntRow,0)=r1;
			rowar(cntRow,1)=r2;
			cntRow++;
			for (j=0;j<2*size;j++)
			{
				tempar(j,0)=c(i,j);
				c(i,j)=c(r2,j);
				c(r2,j)=tempar(j,0);
			}
		}
	}
	double coef ;
	for (i=0;i<size;i++)
	{
		if(i!=(size-1))
		{
			for(j=i+1;j<size;j++)
			{
				coef=c(j,i)/c(i,i);
				for(k=0;k<2*size;k++)
				{
					c(j,k)=c(j,k)-coef*c(i,k);
				}
			}
		}
	}
	for(i=size-1;i>=0;i--)
	{
		for(j=i-1;j>=0;j--)
		{
			coef=c(j,i)/c(i,i);
			for(k=0;k<2*size;k++)
			{
				c(j,k)-=coef*c(i,k);
			}
		}
	}
	for (i=0;i<size;i++)
	{
		coef=c(i,i);
		for(j=0;j<2*size;j++)
		{
			c(i,j)=c(i,j)/coef;
		}
	}
	Matrix<T> b(size,size);
	for(i=0;i<size;i++)
	{
		for(j=size;j<2*size;j++)
		{
			m_Data[i][j-size]=(T)c(i,j);
		}
	}
    
	return true;
    
}
//Invert the matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
Matrix<T>  Matrix<T>::GetInversed() const
{
	Matrix<T> temp(*this);
	temp.Inverse();
	return (temp);
}
//Get the Inverted Matrix
////////////////////////////////////////////////////////////////////////////////////////////


template<class T>
void Matrix<T>::Set_Size(int RowSize, int ColumnSize)
{
	m_RowSize=RowSize;
	m_ColumnSize=ColumnSize;
	m_Data=new T* [m_RowSize];
	for(int i=0;i<m_RowSize;i++)
		m_Data[i]=new T [m_ColumnSize];
    
}
//Set size of the matrix
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
const int & Matrix<T>::Get_RowSize() const
{
	return m_RowSize;
}
//Get Number of rows
////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
const int & Matrix<T>::Get_ColumnSize() const
{
	return m_ColumnSize;
}
//Get Number of columns
////////////////////////////////////////////////////////////////////////////////////////////
template<class T>
Matrix<T> Matrix<T>::operator-(const T& number)
{
	Matrix<T> Temp(*this);
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<m_ColumnSize;j++)
			Temp.Set_Element(i,j,Temp.Get_Element(i,j)-number);
	return Temp;
    
}
//Subtract a number to the matrix
////////////////////////////////////////////////////////////////////////////////////////////
template<class T>
const Matrix<T> & Matrix<T>::operator =(const T&number)
{
	for(int i=0;i<m_RowSize;i++)
		for(int j=0;j<m_ColumnSize;j++)
			m_Data[i][j]=number;
	return *this;
}
//Equalize 2 matrixes
////////////////////////////////////////////////////////////////////////////////////////////
template<class T>
void Matrix<T>::print() const{
	for(int k=0; k<m_RowSize; k++){
		for(int j = 0; j<m_ColumnSize; j++){
			cout << Get_Element(k,j) << ",";
		}
		cout << endl;
	}
}
//Print a Matrix
/////////////////////////////////////////////////////////////////////////////////////////
template<class T>
Matrix<T> Matrix<T>::getKthRow(int k){
	Matrix<T> retMat(1,this->Get_ColumnSize());
    
	if(k < this->Get_RowSize()){
		for(int col = 0; col < this->Get_ColumnSize(); col++){
			retMat.Set_Element(0,col,this->Get_Element(k,col));
		}
	}
	else{
		// ERROR
		cout << "ERROR : k greater than row size of matrix" << endl;
	}
    
	return retMat;
}
// returns k-th row of a given matrix.
/////////////////////////////////////////////////////////////////////////////////////////
template<class T>
void Matrix<T>::setKthRow(int k, Matrix<T> & rowMatrix){
	if(this->Get_ColumnSize() != rowMatrix.Get_ColumnSize() ){
		// ERROR
		cout << "Column sizes of matrices do not match" << endl;
	}
	else if(k > this->Get_RowSize()){
		// ERROR
		cout << "ERROR : k greater than row size of matrix" << endl;
	}
	else{
		for(int col = 0; col < this->Get_ColumnSize(); col++){
			this->Set_Element(k,col,rowMatrix.Get_Element(0,col));
		}
	}
}
// set k-th row of a given matrix.
template<class T>
double Matrix<T>::getRowNorm2(int k){
	double norm = 0.0;
    
	if(k >= this->Get_RowSize()){
		// ERROR
		cout << "k greater than row size of matrix" << endl;
	}
	else{
		for(int i = 0; i < this->Get_ColumnSize(); i++)
			norm = norm + (this->Get_Element(k,i));
	}
    
	return norm;
}
#endif
