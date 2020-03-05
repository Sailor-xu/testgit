//
//  TAMatrix.cpp
//  
//
//  Created by Hongren Chen on 27/09/2017.
//

#include <cstdlib>
#include "TAMatrix.h"

using namespace std;

Matrix::Matrix() {
    rows = 0;
    columns = 0;
    size = 0;
    p = NULL;
}

Matrix::Matrix(int r, int c, double *pt):rows(r), columns(c) {
    size = r * c;
    if (size != 0 && pt != NULL) {
        p = new double[size];
        memcpy(p, pt, size*sizeof(double));
    } else {
        cout<<"unable to create a matrix"<<endl;
    }
}

Matrix::Matrix(int r, int c) :rows(r), columns(c)
{
    size = r*c;
    if (size>0)
    {
        p = new double[size];
        for (unsigned j = 0; j<size; j++) //init
        {
            p[j] = 0.0;
        }
    }
    else
        p = NULL;
}

Matrix::Matrix(int r, int c, double val ) :rows(r), columns(c)// 赋初值val
{
    size = r*c;
    if (size>0)
    {
        p = new double[size];
        for (unsigned j = 0; j<size; j++) //init
        {
            p[j] = val;
        }
    }
    else
        p = NULL;
}

Matrix::Matrix(const Matrix& m) {
    rows = m.rows;
    columns = m.columns;
    size = m.size;
    
    if (rows*columns != 0) {
        p = new double[rows*columns];
        memcpy(p, m.p, rows*columns*sizeof(double));
    } else {
        p = NULL;
    }
}

Matrix::~Matrix() {
    rows = 0;
    columns = 0;
    size = 0;
    if (p != NULL) {
        delete []p;
    }
}

const Matrix& Matrix::operator+() const {
    return *this;
}

const Matrix Matrix::operator-() const {
    double temp[rows*columns];
    for (int i = 0; i < rows*columns; i++) {
        temp[i] = -p[i];
    }
    
    return Matrix(rows, columns, temp);
}

Matrix& Matrix::operator++() {
    for (int i = 0; i < rows*columns; i++) {
        p[i] += 1.0;
    }
    
    return *this;
}

const Matrix Matrix::operator++(int) {
    Matrix old(rows, columns, p);
    for (int i = 0; i < rows*columns; i++) {
        p[i] += 1.0;
    }
    
    return old;
}

Matrix& Matrix::operator--() {
    for (int i = 0; i < rows*columns; i++) {
        p[i] -= 1.0;
    }
    
    return *this;
}

const Matrix Matrix::operator--(int) {
    Matrix old(rows, columns, p);
    for (int i = 0; i < rows*columns; i++) {
        p[i] -= 1.0;
    }
    
    return old;
}

const Matrix Matrix::operator+(const Matrix& m) const {
    if (rows != m.rows || columns != m.columns) {
        cout<<"the sizes of two matrices for operator'+' are different"<<endl;
        exit(0);
    }
    double temp[rows*columns];
    for (int i = 0; i < rows*columns; i++) {
        temp[i] = p[i] + m.p[i];
    }
    
    return Matrix(rows, columns, temp);
}

const Matrix Matrix::operator-(const Matrix& m) const {
    if (rows != m.rows || columns != m.columns) {
        cout<<"the sizes of two matrices for operator'-' are different"<<endl;
        exit(0);
    }
    double temp[rows*columns];
    for (int i = 0; i < rows*columns; i++) {
        temp[i] = p[i] - m.p[i];
    }
    
    return Matrix(rows, columns, temp);
}

const Matrix Matrix::operator*(const Matrix& m) const {
    if (columns != m.rows) {
        cout<<"the sizes of two matrices for operator'*' are not suitable"<<endl;
        exit(0);
    }
    double temp[rows*m.columns];
    memset(temp, 0, rows*m.columns*sizeof(double));
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < m.columns; j++) {
            for (int k = 0; k < columns; k++) {
                temp[i*m.columns+j] += p[i*columns+k]*m.p[j+k*m.columns];
            }
        }
    }
    
    return Matrix(rows, m.columns, temp);
}

const Matrix Matrix::operator+(const double& d) const {
    double temp[rows*columns];
    for (int i = 0; i < rows*columns; i++) {
        temp[i] = p[i] + d;
    }
    
    return Matrix(rows, columns, temp);
}

const Matrix Matrix::operator-(const double& d) const {
    double temp[rows*columns];
    for (int i = 0; i < rows*columns; i++) {
        temp[i] = p[i] - d;
    }
    
    return Matrix(rows, columns, temp);
}

const Matrix Matrix::operator*(const double& d) const {
    double temp[rows*columns];
    for (int i = 0; i < rows*columns; i++) {
        temp[i] = p[i] * d;
    }
    
    return Matrix(rows, columns, temp);
}

const Matrix Matrix::operator/(const double& d) const {
    if (d == 0) {
        cout<<"the divide number is zero"<<endl;
        exit(0);
    }
    double temp[rows*columns];
    for (int i = 0; i < rows*columns; i++) {
        temp[i] = p[i] / d;
    }
    
    return Matrix(rows, columns, temp);
}

const Matrix operator+(const double& d, const Matrix& m) {
    double temp[m.rows*m.columns];
    for (int i = 0; i < m.rows*m.columns; i++) {
        temp[i] = d + m.p[i];
    }
    
    return Matrix(m.rows, m.columns, temp);
}

const Matrix operator-(const double& d, const Matrix& m) {
    double temp[m.rows*m.columns];
    for (int i = 0; i < m.rows*m.columns; i++) {
        temp[i] = d - m.p[i];
    }
    
    return Matrix(m.rows, m.columns, temp);
}

const Matrix operator*(const double& d, const Matrix& m) {
    double temp[m.rows*m.columns];
    for (int i = 0; i < m.rows*m.columns; i++) {
        temp[i] = d * m.p[i];
    }
    
    return Matrix(m.rows, m.columns, temp);
}

Matrix& Matrix::operator=(const Matrix& m) {
    if (this == &m) {
        return *this;
    }
    rows = m.rows;
    columns = m.columns;
    if (p != NULL) {
        delete []p;
    }
    p = new double[m.rows*m.columns];
    memcpy(p, m.p, m.rows*m.columns*sizeof(double));
    
    return *this;
}

Matrix& Matrix::operator+=(const Matrix& m) {
    if (m.p == NULL) {
        return *this;
    }
    if (rows != m.rows || columns != m.columns) {
        cout<<"the sizes of two matrices for operator'+=' are different"<<endl;
        exit(0);
    }
    for (int i = 0; i < rows*columns; i++) {
        p[i] += m.p[i];
    }
    
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& m) {
    if (m.p == NULL) {
        return *this;
    }
    if (rows != m.rows || columns != m.columns) {
        cout<<"the sizes of two matrices for operator'-=' are different"<<endl;
        exit(0);
    }
    for (int i = 0; i < rows*columns; i++) {
        p[i] -= m.p[i];
    }
    
    return *this;
}

Matrix& Matrix::operator*=(const Matrix& m) {
    if (columns != m.rows) {
        cout<<"the sizes of two matrices for operator'*' are not suitable"<<endl;
        exit(0);
    }
    Matrix old(rows, columns, p);
    delete []p;
    p = new double[rows*m.columns];
    memset(p, 0, rows*m.columns*sizeof(double));
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < m.columns; j++) {
            for (int k = 0; k < columns; k++) {
                p[i*m.columns+j] += old.p[i*columns+k]*m.p[j+k*m.columns];
            }
        }
    }
    
    return *this;
}

Matrix& Matrix::operator+=(const double& d) {
    if (d == 0) {
        return *this;
    }
    for (int i = 0; i < rows*columns; i++) {
        p[i] += d;
    }
    
    return *this;
}

Matrix& Matrix::operator-=(const double& d) {
    if (d == 0) {
        return *this;
    }
    for (int i = 0; i < rows*columns; i++) {
        p[i] -= d;
    }
    
    return *this;
}

Matrix& Matrix::operator*=(const double& d) {
    for (int i = 0; i < rows*columns; i++) {
        p[i] *= d;
    }
    
    return *this;
}

Matrix& Matrix::operator/=(const double& d) {
    if (d == 0) {
        cout<<"the divide number is zero"<<endl;
        exit(0);
    }
    for (int i = 0; i < rows*columns; i++) {
        p[i] /= d;
    }
    
    return *this;
}

bool Matrix::operator==(const Matrix& m) const {
    if (rows != m.rows || columns != m.columns) {
        return false;
    }
    for (int i = 0; i < rows*columns; i++) {
        if (p[i] != m.p[i]) {
            return false;
        }
    }
    return true;
}

ostream& operator<<(ostream& os, const Matrix& m) {
    for (int i = 0; i < m.rows; i++) {
        for (int j = 0; j < m.columns; j++) {
            os<<m.p[i*m.columns+j]<<' ';
        }
        os<<endl;
    }
    
    return os;
}

double& Matrix::takeElem(int row, int column) const {
    return *(p + (row-1) * columns + (column-1));
}

Matrix Matrix::takeRow(int row) const {
    double temp[columns];
    for (int i = 0; i < columns; i++) {
        temp[i] = p[(row-1)*columns+i];
    }
    
    return Matrix(1, columns, temp);
}

Matrix Matrix::takeColumn(int column) const {
    double temp[rows];
    for (int i = 0; i < rows; i++) {
        temp[i] = p[column-1+columns*i];
    }
    
    return Matrix(rows, 1, temp);
}
int& Matrix::getRows() const {
    return rows;
}
int& Matrix::getColumns() const {
    return columns;
}

Matrix Matrix::takeSubMatrix(int row1, int column1, int row2, int column2) const {
    if (row1 < 1 || row1 > rows || row2 < 1 || row2 > rows ||
        column1 < 1 || column1 > columns || column2 < 1 || column2 > columns) {
        cout<<"out of range of matrix"<<endl;
        exit(0);
    }
    if (row1 > row2) {
        int t = row1;
        row1 = row2;
        row2 = t;
    }
    if (column1 > column2) {
        int t = column1;
        column1 = column2;
        column2 = t;
    }
    
    double temp[(row2-row1+1)*(column2-column1+1)];
    for (int i = row1; i <= row2; i++) {
        for (int j = column1; j <= column2; j++) {
            temp[(i-row1)*(column2-column1+1)+(j-column1)] = p[(i-1)*columns+(j-1)];
        }
    }
    
    return Matrix(row2-row1+1, column2-column1+1, temp);
}

Matrix Matrix::inverse() const {
    if(rows < 2 && columns < 2){
        Matrix mat(rows, columns, 1 / p[0]);
        return mat;
    } else {
        //计算伴随矩阵
        double adjoint_matrix[rows * columns];
        double det[(rows - 1) * (columns - 1)];
        for (int m = 0; m < rows; m++) {
            for (int n = 0; n < columns; n++) {
                double k = pow(-1, m + n);
                for (int i = 0; i < rows; i++) {
                    if (i == m) {
                        continue;
                    } else if (i < m) {
                        for (int j = 0; j < columns; j++) {
                            if (j == n) {
                                continue;
                            } else if (j < n) {
                                det[i * (columns - 1) + j] = p[i * columns + j];
                            } else {
                                det[i * (columns - 1) + (j - 1)] = p[i * columns + j];
                            }
                        }
                    } else {
                        for (int j = 0; j < columns; j++) {
                            if (j == n) {
                                continue;
                            } else if (j < n) {
                                det[(i - 1) * (columns - 1) + j] = p[i * columns + j];
                            } else {
                                det[(i - 1) * (columns - 1) + (j - 1)] = p[i * columns + j];
                            }
                        }
                    }
                }
                Matrix detMat(rows - 1, columns - 1, det);
                adjoint_matrix[n * columns + m] = k * detMat.detMatrix();
            }
        }
        Matrix mat(rows, columns, adjoint_matrix);
        //计算逆矩阵
        return mat / this->detMatrix();
    }
}

Matrix Matrix::transpose() const {
    double temp[columns*rows];
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            temp[j*rows+i] = p[i*columns+j];
        }
    }
    
    return Matrix(columns, rows, temp);
}

double Matrix::detMatrix() const {
    if (rows != columns) {
        cout<<"the matrix has no detMatrix"<<endl;
        exit(0);
    }
    if (rows <= 0) {
        cout<<"the matrix size is incorrect"<<endl;
        exit(0);
    }
    else if (rows == 1) {
        return p[0];
    }
    else if (rows == 2) {
        return (p[0]*p[3] - p[1]*p[2]);
    }
    else if (rows == 3) {
        return (p[0]*p[4]*p[8] + p[1]*p[5]*p[6] + p[2]*p[3]*p[7] - p[2]*p[4]*p[6] - p[1]*p[3]*p[8]- p[0]*p[5]*p[7]);
    }
    else {
        double det = 0;  //降阶法
        for (int i = 0; i < columns; i++) {
            double m = p[i] * pow(-1, i);  //计算系数
            double temp[(rows-1)*(columns-1)];
            //计算该系数对应的余子式
            for (int j = 1; j < rows; j++) {
                for (int k = 0; k < columns; k++) {
                    if (k == i) {
                        continue;
                    }
                    else if (k < i) {
                        temp[(j-1)*(columns-1)+k] = p[j*columns+k];
                    }
                    else {
                        temp[(j-1)*(columns-1)+k-1] = p[j*columns+k];
                    }
                }
            }
            Matrix mat(rows-1, columns-1, temp);
            //累加
            det += m * mat.detMatrix();
        }
        return det;
    }
}

Matrix Matrix::eye() {
    for (unsigned i = 0; i< rows; i++) {
        for (unsigned j = 0; j < columns; j++) {
            if (i == j) {
                p[i * columns + j] = 1.0;
            } else {
                p[i * columns + j] = 0.0;
            }
        }
    }
    return *this;
}

Matrix Matrix::diag()
{
    if (rows != columns)
    {
        Matrix m(0,0);
        cout << "diag():row != col" << endl;
        return m;
    }
    Matrix m(rows, rows);
    for (unsigned i = 0; i<rows; i++)
    {
        m.p[i*rows + i] = p[i*rows + i];
    }
    return m;
}

Matrix Matrix::cov(bool flag)
{
    if (columns == 0)
    {
        Matrix m(0, 0);
        return m;
    }
    double *mean = new double[columns]; //均值向量

    for (unsigned j = 0; j<columns; j++) //init
    {
        mean[j] = 0.0;
    }
    Matrix ret(rows, columns);
    for (unsigned j = 0; j<columns; j++) //mean
    {
        for (unsigned i = 0; i<rows; i++)
        {
            mean[j] += p[i*columns + j];
        }
        mean[j] /= rows;
    }
    unsigned i, k, j;
    for (i = 0; i<columns; i++) //第一个变量
    {
        for (j = i; j<columns; j++) //第二个变量
        {
            for (k = 0; k<rows; k++) //计算
            {
                ret.p[i * columns + j] += (p[k*columns + i] - mean[i])*(p[k*columns + j] - mean[j]);

            }
            if (flag)
            {
                ret.p[i * columns + j] /= (rows-1);
            }
            else
            {
                ret.p[i * columns + j] /= (rows);
            }
        }
    }
    for (i = 0; i<columns; i++) //补全对应面
    {
        for (j = 0; j<i; j++)
        {
            ret.p[i * columns + j] = ret.p[j * columns + i];
        }
    }
    return ret;
}

void  Matrix::QR(Matrix &Q, Matrix &R) const
{
    //如果A不是一个二维方阵，则提示错误，函数计算结束
    if (rows != columns)
    {
        printf("ERROE: QR() parameter A is not a square matrix!\n");
        return;
    }
    const unsigned N = unsigned(rows);
    double *a = new double[N];
    double *b = new double[N];

    for (unsigned j = 0; j < N; ++j)  //(Gram-Schmidt) 正交化方法
    {
        for (unsigned i = 0; i < N; ++i)  //第j列的数据存到a，b
            a[i] = b[i] = p[i * N + j];

        for (unsigned i = 0; i<j; ++i)  //第j列之前的列
        {
            R.p[i * N + j] = 0;  //
            for (unsigned m = 0; m < N; ++m)
            {
                R.p[i * N + j] += a[m] * Q.p[m *N + i]; //R[i,j]值为Q第i列与A的j列的内积
            }
            for (unsigned m = 0; m < N; ++m)
            {
                b[m] -= R.p[i * N + j] * Q.p[m * N + i]; //
            }
        }

        double norm = 0;
        for (unsigned i = 0; i < N; ++i)
        {
            norm += b[i] * b[i];
        }
        norm = sqrt(norm);

        R.p[j*N + j] = norm; //向量b[]的2范数存到R[j,j]

        for (unsigned i = 0; i < N; ++i)
        {
            Q.p[i * N + j] = b[i] / norm; //Q 阵的第j列为单位化的b[]
        }
    }
    delete[]a;
    delete[]b;
}

Matrix Matrix::eig_val(unsigned _iters)
{
    if (size == 0 || rows != columns)
    {
        cout << "矩阵为空或者非方阵！" << endl;
        Matrix rets(0,0);
        return rets;
    }
    const unsigned N = unsigned(rows);
    Matrix matcopy(*this);//备份矩阵
    Matrix Q(N,N), R(N,N);
    /*当迭代次数足够多时,A 趋于上三角矩阵，上三角矩阵的对角元就是A的全部特征值。*/
    for (unsigned k = 0; k < _iters; ++k)
    {
        QR(Q, R);
        *this = R*Q;
    }
    Matrix val = diag();
    *this = matcopy;//恢复原始矩阵；
    return val;
}
Matrix Matrix::eig_vect(unsigned _iters)
{
    if (size == 0 || rows != columns)
    {
        cout << "矩阵为空或者非方阵！" << endl;
        Matrix rets(0,0);
        return rets;
    }
    if (detMatrix() == 0)
    {
        cout << "非满秩矩阵没法用QR分解计算特征向量！" << endl;
        Matrix rets(0,0);
        return rets;
    }
    const unsigned N = unsigned(rows);
    Matrix matcopy(*this);//备份矩阵
    Matrix Q(N,N), R(N,N);
    /*当迭代次数足够多时,A 趋于上三角矩阵，上三角矩阵的对角元就是A的全部特征值。*/
    for (unsigned k = 0; k < _iters; ++k)
    {
        QR(Q, R);
    }
    Matrix ret = Q;

    *this = matcopy;//恢复原始矩阵；
    return ret;
}
/*Matrix Matrix::eig_vect(unsigned _iters)
{
    if (size == 0 || rows != columns)
    {
        cout << "矩阵为空或者非方阵！" << endl;
        Matrix rets(0,0);
        return rets;
    }
    if (detMatrix() == 0)
    {
        cout << "非满秩矩阵没法用QR分解计算特征向量！" << endl;
        Matrix rets(0,0);
        return rets;
    }
    Matrix matcopy(*this);//备份矩阵
    Matrix eigenValue = eig_val(_iters);
    Matrix ret(rows, rows);
    const int NUM = columns;
    double eValue;
    double sum, midSum, diag;
    Matrix copym(*this);
    for (unsigned count = 0; count < NUM; ++count)
    {
        //计算特征值为eValue，求解特征向量时的系数矩阵
        *this = copym;
        eValue = eigenValue.p[count * rows + count];

        for (unsigned i = 0; i < columns; ++i)//A-lambda*I
        {
            p[i * columns + i] -= eValue;
        }
        //cout<<*this<<endl;
        //将 this为阶梯型的上三角矩阵
        for (unsigned i = 0; i < rows - 1; ++i)
        {
            diag = p[i*columns + i];  //提取对角元素
            for (unsigned j = i; j < columns; ++j)
            {
                p[i*columns + j] /= diag; //【i,i】元素变为1
            }
            for (unsigned j = i + 1; j<rows; ++j)
            {
                diag = p[j *  columns + i];
                for (unsigned q = i; q < columns; ++q)//消去第i+1行的第i个元素
                {
                    p[j*columns + q] -= diag*p[i*columns + q];
                }
            }
        }
        //特征向量最后一行元素置为1
        midSum = ret.p[(ret.rows - 1) * ret.columns + count] = 1;
        for (int m = rows - 2; m >= 0; --m)
        {
            sum = 0;
            for (unsigned j = m + 1; j < columns; ++j)
            {
                sum += p[m *  columns + j] * ret.p[j * ret.columns + count];
            }
            sum = -sum / p[m *  columns + m];
            midSum += sum * sum;
            ret.p[m * ret.columns + count] = sum;
        }
        midSum = sqrt(midSum);
        for (unsigned i = 0; i < ret.rows; ++i)
        {
            ret.p[i * ret.columns + count] /= midSum; //每次求出一个列向量
        }
    }
    *this = matcopy;//恢复原始矩阵；
    return ret;
}*/



//构造函数，构建空矩阵
CMatrix::CMatrix(){
    row=0;column=0;data=nullptr;
}
//由Complex转化为1x1矩阵
CMatrix::CMatrix(std::complex item){
    row=1;column=1;this->data=new std::complex[row*column];
    data[0]=item;
}
//由实矩阵转化为复数矩阵
CMatrix::CMatrix(Matrix& a){
    row=a.getRows();
    column=a.getColumns();
    this->data=new std::complex[row*column];
    for (int i=0;i<row;i++)
        for(int j=0;j<column;j++)
            data[i*column+j]=std::complex(a.takeElem(i, j),0);
}
//由实部和虚部组成复数矩阵
CMatrix::CMatrix(Matrix& real,Matrix& imag){
    row=real.getRows();
    column=real.getColumns();
    this->data=new std::complex[row*column];
    for (int i=0;i<row;i++)
        for(int j=0;j<column;j++)
        data[i*column+j]=std::complex(real.takeElem(i, j),imag.takeElem(i, j));
}

//构建矩阵，不赋值
CMatrix::CMatrix(int row,int column){
    this->row=row;
    this->column=column;
    this->data=new std::complex[row*column];
}
//构建值矩阵，矩阵中每个元素都等于给定值
CMatrix::CMatrix(int row,int column,std::complex data){
    this->row=row;
    this->column=column;
    this->data=new std::complex[row*column];
    for (int i=0;i<row*column;i++)
        this->data[i]=data;
}
//根据数组构建矩阵,复制元素
CMatrix::CMatrix(int row,int column,double *x,double *y){
    this->row=row;
    this->column=column;
    this->data=new std::complex[row*column];
    for (int i=0;i<row*column;i++)
        this->data[i]=std::complex(x[i],y[i]);
}
//根据数组构建矩阵,复制元素
CMatrix::CMatrix(double *x,double *y,int row,int column){
    this->row=row;
    this->column=column;
    this->data=new std::complex[row*column];
    for (int i=0;i<row*column;i++)
        this->data[i]=std::complex(x[i],y[i]);
}
//根据数组构建矩阵,复制元素
CMatrix::CMatrix(int row,int column,std::complex *data){
    this->row=row;
    this->column=column;
    this->data=new std::complex[row*column];
    for (int i=0;i<row*column;i++)
        this->data[i]=data[i];
}
//根据数组构建矩阵,复制元素
CMatrix::CMatrix(std::complex *data,int row,int column){
    this->row=row;
    this->column=column;
    this->data=new std::complex[row*column];
    for (int i=0;i<row*column;i++)
        this->data[i]=data[i];
}
//根据向量构建矩阵
CMatrix::CMatrix(std::vector<std::complex> a){
    this->row=int(a.size());
    this->column=1;
    this->data=new std::complex[row*column];
    for (int i=0;i<a.size();i++)
        this->data[i] = a[i];
}
//根据向量构建矩阵
CMatrix::CMatrix(std::vector<std::complex> a,int row,int column){
    if (a.size()<row*column)
    {
        std::cout<<"无法根据向量vector生成矩阵"<<std::endl;
    return;
    //*this=NULL;
    }
    this->row=row;
    this->column=column;
    this->data=new std::complex[row*column];
    for (int i=0;i<a.size();i++)
    this->data[i] = a[i];
}
//根据向量构建矩阵
CMatrix::CMatrix(std::vector<std::vector<std::complex>> a){
    if (a.size()>0 && a[0].size()>0)
    {
        this->row=int(a.size());
        this->column=a[0].size();
        this->data=new std::complex[row*column];
        for (int i=0;i<row;i++)
            for (int j=0;j<column;j++)
                this->data[i*column+j] = a[i][j];
    }
}
//复制构造函数
CMatrix::CMatrix(const CMatrix& another){
    this->row=another.row;
    this->column=another.column;
    this->data=new std::complex[row*column];
    for (int i=0;i<row*column;i++)
        this->data[i]=another.data[i];
}
//重载赋值运算符，保留在matlab中赋值等于深复制的效果
CMatrix& CMatrix::operator=(const CMatrix& another){
    if(this==&another)
        return *this;
    this->row=another.row;
    this->column=another.column;
    if(data !=nullptr)   //必须新释放之前的内存，不然内存泄漏了
        delete[] data;
    this->data=new std::complex[row*column];
    for (int i=0;i<row*column;i++)
        this->data[i]=another.data[i];
    return *this;
}
//重载()运算符,下标从0开始的习惯，读取时调用
std::complex CMatrix::operator()(int i,int j) const
{
    if (i>this->row-1 || i<0)
        std::cout<<"读取矩阵元素，超出矩阵范围"<<std::endl;
    if (j>this->column-1 || j<0)
        std::cout<<"读取矩阵元素，超出矩阵范围"<<std::endl;
    return data[i*column+j];
}
//重载()运算符,下标从0开始的习惯，写入时调用
std::complex& CMatrix::operator()(int i,int j)
{
    if (i>this->row-1 || i<0)
        std::cout<<"读取矩阵元素，超出矩阵范围"<<std::endl;
    if (j>this->column-1 || j<0)
        std::cout<<"读取矩阵元素，超出矩阵范围"<<std::endl;
    return data[i*column+j];
}
//重载()运算符,下标从0开始的习惯，读取时调用
std::complex CMatrix::operator()(int i) const
{
    if (i>this->row*this->column-1 || i<0)
        std::cout<<"读取矩阵元素，超出矩阵范围"<<std::endl;
    return data[i];
}
//重载()运算符,下标从0开始的习惯
std::complex& CMatrix::operator()(int i){
    if (i>this->row*this->column-1 || i<0)
        std::cout<<"读取矩阵元素，超出矩阵范围"<<std::endl;
    return data[i];
}

//重载加法运算符
CMatrix CMatrix::operator+(std::complex xi){
    CMatrix p(row,column);
    for (int i=0;i<p.row*p.column;i++)
        p.data[i] = data[i]+xi;
    return p;
}
//重载加法运算符
CMatrix& CMatrix::operator+=(std::complex xi){
    for (int i=0;i<row*column;i++)
        this->data[i] += xi;
    return *this;
}
//重载加法运算符
CMatrix CMatrix::operator+(CMatrix& another){
    if(!another.data)
        return *this;
    if (row!=another.row || column!=another.column)
        std::cout<<"矩阵不匹配"<<std::endl;
    CMatrix p(row,column);
    for (int i=0;i<row*column;i++)
        p(i)=data[i]+another(i);
    return p;
}
//重载加法运算符
CMatrix& CMatrix::operator+=(CMatrix& another){
    if(!another.data)
        return *this;
    if (row!=another.row || column!=another.column)
        std::cout<<"矩阵不匹配"<<std::endl;
    for (int i=0;i<row*column;i++)
        this->data[i]+=another.data[i];
    return *this;
}
//重载减法运算符
CMatrix CMatrix::operator-(std::complex xi){
    CMatrix p(row,column);
    for (int i=0;i<p.row*p.column;i++)
        p.data[i] = data[i]-xi;
    return p;
}
//重载减法运算符
CMatrix& CMatrix::operator-=(std::complex xi){
    for (int i=0;i<row*column;i++)
        this->data[i] -= xi;
    return *this;
}
//重载减法运算符
CMatrix CMatrix::operator-(CMatrix& another){
    if(!another.data)
        return *this;
    if (row!=another.row || column!=another.column)
        std::cout<<"矩阵不匹配"<<std::endl;
    CMatrix p(row,column);
    for (int i=0;i<row*column;i++)
        p.data[i]=data[i]-another.data[i];
    return p;
}
//重载减法运算符
CMatrix& CMatrix::operator-=(CMatrix& another){
    if(!another.data)
        return *this;
    if (row!=another.row || column!=another.column)
        std::cout<<"矩阵不匹配"<<std::endl;
    for (int i=0;i<row*column;i++)
        this->data[i]-=another.data[i];
    return *this;
}
//重载乘法运算符
CMatrix CMatrix::operator*(std::complex xi){
    CMatrix p(row,column);
    for (int i=0;i<p.row*p.column;i++)
        p.data[i] = this->data[i]*xi;
    return p;
}
//重载乘法运算符
CMatrix& CMatrix::operator*=(std::complex xi){
    for (int i=0;i<row*column;i++)
        this->data[i] *= xi;
    return *this;
}
//重载乘法运算符
CMatrix CMatrix::operator*(CMatrix& another){
    if(!another.data || !this->data)
        std::cout<<"矩阵乘法数据为空"<<std::endl;
    if (column!=another.row)
        std::cout<<"矩阵不匹配"<<std::endl;
    CMatrix p(row,another.column);
    for (int i=0;i<p.row;i++)
        for (int j=0;j<p.column;j++)
        {
            std::complex sum=0;
            for (int t=0;t<column;t++)  //第一个矩阵的第i行乘以第二个矩阵的第j列
                sum+=data[i*column+t]*another(t,j);
            p(i,j) = sum;
        }
    return p;
}
//重载除法运算符
CMatrix CMatrix::operator/(std::complex xi){
    CMatrix p(row,column);
    for (int i=0;i<p.row*p.column;i++)
        p.data[i] = data[i]/xi;
    return p;
}
//重载除法运算符
CMatrix& CMatrix::operator/=(std::complex xi){
    for (int i=0;i<row*column;i++)
        this->data[i] /= xi;
    return *this;
}
//重载矩阵等号==运算符
bool CMatrix::operator==(CMatrix& another){
    if (this->row!= another.row)
        return false;
    if (this->column != another.column)
        return false;

    for (int i = 0; i < row*column; i ++)
        if (FloatNotEqual(this->data[i].real(),another.data[i].real()) || FloatNotEqual(this->data[i].imag(),another.data[i].imag()))
            return false;
    return true;
}
//重载矩阵等号!=运算符
bool CMatrix::operator!=(CMatrix& another){
    if (this->row!= another.row)
        return true;
    if (this->column != another.column)
        return true;

    for (int i = 0; i < row*column; i ++)
        if (FloatNotEqual(this->data[i].real(),another.data[i].real()) || FloatNotEqual(this->data[i].imag(),another.data[i].imag()))
            return true;
    return false;
}
//重载一元操作符：正号
CMatrix CMatrix::operator +() const
{
    return *this;           //不用作任何处理，维持原状
}
//重载一元操作符：负号
CMatrix CMatrix::operator -() const
{
    CMatrix p(row,column);
    for (int i=0;i<row*column;i++)
        p.data[i]=-this->data[i];
    return p;
}
//矩阵转置，不改变源矩阵
CMatrix CMatrix::transpose()
{
    CMatrix p(column,row);
    for(int i=0;i<row;i++)
        for(int j=0;j<column;j++)
            p(j,i)=(*this)(i,j);
    return p;
}

//矩阵左扩展，改变了源矩阵
void CMatrix::append_left(CMatrix& a){
    if (a.row==0)
        return;
    if (this->row==0){
        this->row = a.row;
        this->column=a.column;
        this->data = new std::complex[row*column];
        for (int i=0;i<row*column;i++)
            this->data[i]=a.data[i];
        return;
    }
    int new_column=a.column+column;
    if (a.row==this->row)
    {
        std::complex* datanew=new std::complex[row*new_column];
        for (int i=0;i<a.row;i++)
        {
            for (int j=0;j<a.column;j++)
                datanew[i*new_column+j]=a(i,j);
            for (int j=0;j<column;j++)
                datanew[i*new_column+j+a.column]=this->data[i*column+j];
        }
        this->column=new_column;
        delete[] this->data;
        this->data = datanew;
    }else
    {
        std::cout<<"合并矩阵行数不等"<<std::endl;
    }
}
//向量首部添加元素
void CMatrix::append_left(std::complex x){
    CMatrix a(1,1,x);
    if(column==1)
        this->append_top(a);
    else if(row==1)
        this->append_left(a);
}
//矩阵右扩展，改变了源矩阵
void CMatrix::append_right(CMatrix& a){
    if (a.row==0)
        return;
    if (this->row==0){
        this->row = a.row;
        this->column=a.column;
        this->data = new std::complex[row*column];
        for (int i=0;i<row*column;i++)
            this->data[i]=a.data[i];
        return;
    }
    int new_column=a.column+column;
    if (a.row==this->row)
    {
        std::complex* datanew=new std::complex[row*new_column];
        for (int i=0;i<row;i++)
        {
            for (int j=0;j<column;j++)
                datanew[i*new_column+j]=this->data[i*column+j];
            for (int j=0;j<a.column;j++)
                datanew[i*new_column+j+column]=a(i,j);
        }
        this->column=new_column;
        delete[] this->data;
        this->data = datanew;
    }else
    {
        std::cout<<"合并矩阵行数不等"<<std::endl;
    }
}
//向量尾部添加元素
void CMatrix::append_right(std::complex x){
    CMatrix a(1,1,x);
    if(column==1)
        this->append_bottom(a);
    else if(row==1)
        this->append_right(a);
}
//矩阵上扩展，改变了源矩阵
void CMatrix::append_top(CMatrix& a){
    if (a.row==0)
        return;
    if (this->row==0){
        this->row = a.row;
        this->column=a.column;
        this->data = new std::complex[row*column];
        for (int i=0;i<row*column;i++)
            this->data[i]=a.data[i];
        return;
    }
    int new_row=a.row+row;
    if (a.column==this->column)
    {
        std::complex* datanew=new std::complex[new_row*column];
        for (int i=0;i<a.row;i++)
            for (int j=0;j<a.column;j++)
                datanew[i*column+j]=a(i,j);
        for (int i=0;i<row;i++)
            for (int j=0;j<column;j++)
                datanew[(i+a.row)*column+j]=this->data[i*column+j];
        this->row=new_row;
        delete[] this->data;
        this->data = datanew;
    }else
    {
        std::cout<<"合并矩阵列数不等"<<std::endl;
    }
}
//向量首部添加元素
void CMatrix::append_top(std::complex x){
    CMatrix a(1,1,x);
    if(column==1)
        this->append_top(a);
    else if(row==1)
        this->append_left(a);
}
//矩阵下扩展，改变了源矩阵
void CMatrix::append_bottom(CMatrix& a){
    if (a.row==0)
        return;
    if (this->row==0){
        this->row = a.row;
        this->column=a.column;
        this->data = new std::complex[row*column];
        for (int i=0;i<row*column;i++)
            this->data[i]=a.data[i];
        return;
    }
    int new_row=a.row+row;
    if (a.column==this->column)
    {
        std::complex* datanew=new std::complex[new_row*column];
        for (int i=0;i<row;i++)
            for (int j=0;j<column;j++)
                datanew[i*column+j]=this->data[i*column+j];
        for (int i=0;i<a.row;i++)
            for (int j=0;j<a.column;j++)
                datanew[(i+row)*column+j]=a(i,j);
        this->row=new_row;
        delete[] this->data;
        this->data = datanew;
    }else
    {
        std::cout<<"合并矩阵列数不等"<<std::endl;
    }
}
//向量尾部添加元素
void CMatrix::append_bottom(std::complex x){
    CMatrix a(1,1,x);
    if(column==1)
        this->append_bottom(a);
    else if(row==1)
        this->append_right(a);
}
//删除指定行，改变了源矩阵
void CMatrix::remove_row(int rowindex)  //删除指定行
{
    if (rowindex>=row || rowindex<0)
    {
        std::cout<<"删除行范围越界"<<std::endl;
        return;
    }
    std::complex* datetemp = data;
    int newrow=row-1;
    data = new std::complex[newrow*column];
    for (int i=0;i<rowindex;i++)
        for (int j=0;j<this->column;j++)
            this->data[i*column+j] = datetemp[i*column+j];
    for (int i=rowindex;i<newrow;i++)
        for (int j=0;j<this->column;j++)
            this->data[i*column+j] = datetemp[(i+1)*column+j];
    row=newrow;
    delete[] datetemp;
}
//删除指定行，改变了源矩阵
void CMatrix::remove_row(int row1,int row2)  //删除多行,包含两个下标
{
    if (row1>=row || row2>=row || row1<0||row2<0)
    {
        std::cout<<"删除行范围越界"<<std::endl;
        return;
    }
    std::complex* datetemp = data;
    int newrow=row-row2+row1-1;
    data = new std::complex[newrow*column];
    for (int i=0;i<row1;i++)
        for (int j=0;j<this->column;j++)
            this->data[i*column+j] = datetemp[i*column+j];
    for (int i=row1;i<newrow;i++)
        for (int j=0;j<this->column;j++)
            this->data[i*column+j] = datetemp[(i+row2-row1+1)*column+j];
    row=newrow;
    delete[] datetemp;
}
//删除指定列，改变了源矩阵
void CMatrix::remove_column(int columnindex)  //删除指定列
{
    if (columnindex>=column || columnindex<0)
    {
        std::cout<<"删除列范围越界"<<std::endl;
        return;
    }
    std::complex* datetemp = data;
    int newcolumn=column-1;
    data = new std::complex[row*newcolumn];
    for (int i=0;i<row;i++)
    {
        for (int j=0;j<columnindex;j++)
            this->data[i*newcolumn+j] = datetemp[i*column+j];
        for (int j=columnindex;j<newcolumn;j++)
            this->data[i*newcolumn+j] = datetemp[i*column+j+1];
    }
    column=newcolumn;
    delete[] datetemp;
}
//删除指定列，改变了源矩阵
void CMatrix::remove_column(int column1,int column2)  //删除多列,包含两个下标
{
    if (column1>=column || column2>=column || column1<0||column2<0)
    {
        std::cout<<"删除列范围越界"<<std::endl;
        return;
    }
    std::complex* datetemp = data;
    int newcolumn=column-column2+column1-1;
    data = new std::complex[row*newcolumn];

    for (int i=0;i<row;i++)
    {
        for (int j=0;j<column1;j++)
            this->data[i*newcolumn+j] = datetemp[i*column+j];
        for (int j=column1;j<newcolumn;j++)
            this->data[i*newcolumn+j] = datetemp[i*column+j+column2-column1+1];
    }
    column=newcolumn;
    delete[] datetemp;
}
//删除指定的行和列
void CMatrix::remove_row_column(int rowindex,int columnindex)
{
    if (columnindex>=column || columnindex<0 || rowindex<0 || rowindex>=row)
    {
        std::cout<<"删除行列范围越界"<<std::endl;
        return;
    }
    std::complex* datetemp = data;
    int newcolumn=column-1;
    int newrow=row-1;
    data = new std::complex[newrow*newcolumn];
    for (int i=0;i<rowindex;i++)
    {
        for (int j=0;j<columnindex;j++)
            this->data[i*newcolumn+j] = datetemp[i*column+j];
        for (int j=columnindex;j<newcolumn;j++)
            this->data[i*newcolumn+j] = datetemp[i*column+j+1];
    }
    for (int i=rowindex;i<newrow;i++)
    {
        for (int j=0;j<columnindex;j++)
            this->data[i*newcolumn+j] = datetemp[(i+1)*column+j];
        for (int j=columnindex;j<newcolumn;j++)
            this->data[i*newcolumn+j] = datetemp[(i+1)*column+j+1];
    }
    column=newcolumn;
    row =newrow;
    delete[] datetemp;
}
//删除多行和多列
void CMatrix::remove_row_column(int row1,int row2,int column1,int column2)
{
    if (column1>=column || column1<0 || row1<0 || row1>=row || column2>=column || column2<0 || row2<0 || row2>=row)
    {
        std::cout<<"删除行列范围越界"<<std::endl;
        return;
    }
    std::complex* datetemp = data;
    int newcolumn=column-column2+column1-1;
    int newrow=row-row2+row1-1;
    data = new std::complex[newrow*newcolumn];
    for (int i=0;i<row1;i++)
    {
        for (int j=0;j<column1;j++)
            this->data[i*newcolumn+j] = datetemp[i*column+j];
        for (int j=column1;j<newcolumn;j++)
            this->data[i*newcolumn+j] = datetemp[i*column+j+column2-column1+1];
    }
    for (int i=row1;i<newrow;i++)
    {
        for (int j=0;j<column1;j++)
            this->data[i*newcolumn+j] = datetemp[(i+row2-row1+1)*column+j];
        for (int j=column1;j<newcolumn;j++)
            this->data[i*newcolumn+j] = datetemp[(i+row2-row1+1)*column+j+column2-column1+1];
    }
    column=newcolumn;
    row =newrow;
    delete[] datetemp;
}
//替换局部矩阵，改变了源矩阵，
void CMatrix::replace(CMatrix &a,int row1,int column1)
{
    int endrow=std::min(this->row,a.row+row1);
    int endcolumn=std::min(this->column,a.column+column1);
    for (int i=row1;i<endrow;i++)
        for (int j=column1;j<endcolumn;j++)
            this->data[i*column+j] = a((i-row1),j-column1);
}

//点运算
CMatrix CMatrix::dot(CMatrix& a,std::string operation)
{
    if (row!=a.row || column!=a.column)
    {
        std::cout<<"点运算矩阵不匹配"<<std::endl;
        return NULL;
    }
    CMatrix p(row,column);
    if (operation=="+")
        for (int i=0;i<row;i++)
            for(int j=0;j<column;j++)
                p(i,j)=(*this)(i,j)+a(i,j);
    else if (operation=="-")
        for (int i=0;i<row;i++)
            for(int j=0;j<column;j++)
                p(i,j)=(*this)(i,j)-a(i,j);
    else if (operation=="*")
        for (int i=0;i<row;i++)
            for(int j=0;j<column;j++)
                p(i,j)=(*this)(i,j)*a(i,j);
    else if (operation=="/")
        for (int i=0;i<row;i++)
            for(int j=0;j<column;j++)
                p(i,j)=(*this)(i,j)/a(i,j);
    else if (operation=="\\")
        for (int i=0;i<row;i++)
            for(int j=0;j<column;j++)
                p(i,j)=a(i,j)/(*this)(i,j);
    return p;
}

//查询是否是向量
bool CMatrix::isVector()
{
    if (row==1||column==1)
        return true;
    return false;
}
//实部矩阵
Matrix CMatrix::real(){
    double *p = new double[row * column];
    for (int i=0;i<row * column;i++){
        p[i] = data[i].real();
    }
    Matrix ret(row, column, p);
    return ret;
}
//虚部矩阵
Matrix CMatrix::imag(){
    double *p = new double[row * column];
    for (int i=0;i<row * column;i++){
        p[i] = data[i].imag();
    }
    Matrix ret(row, column, p);
    return ret;
}
//模值矩阵
Matrix CMatrix::abs(){
    double *p = new double[row * column];
    for (int i=0;i<row * column;i++)
        p[i]=(double)std::sqrt((this->data[i].imag())*(this->data[i].imag())+(this->data[i].real())*(this->data[i].real()));
    Matrix ret(row, column, p);
    return ret;
}
//角度矩阵
Matrix CMatrix::angle(){
    double *p = new double[row * column];
    for (int i=0;i<row * column;i++)
        p[i]=(double)std::atan((this->data[i].imag())/(this->data[i].real()));
    Matrix ret(row, column, p);
    return ret;
}
//共轭矩阵
CMatrix CMatrix::conj(){
    CMatrix p(row,column);
    for (int i=0;i<p.row*p.column;i++)
        p(i)=std::complex(this->data[i].real(),-(this->data[i].imag()));
    return p;
}

//按行优先转化为series
std::complex* CMatrix::toSeries(){
    std::complex* back = new std::complex[row*column];
    for (int i=0;i<column*row;i++)
        back[i] = data[i];
    return back;
}
//转化为二位数组
std::complex** CMatrix::toArray(){
    std::complex** back = new std::complex*[column];
    int n=0;
    for (int i=0;i<row;i++){
        back[i]=new std::complex[column];
        for(int j=0;j<column;j++)
            back[i][j] = data[i*column+j];
    }
    return back;
}
//按行优先转化为一维vector
std::vector<std::complex> CMatrix::toVector(){
    int t=0;
    std::vector<std::complex> back(row*column);
    for (int i=0;i<row;i++)
        for (int j=0;j<column;j++)
            back[t++]=data[i*column+j];
    return back;
}
//按行优先转化为二维vector
std::vector<std::vector<std::complex>> CMatrix::toVector2(){
    int rowindex=0;
    int columnindex=0;
    std::vector<std::vector<std::complex>> back(row);
    for (int i=0;i<row;i++)
    {
        columnindex=0;
        std::vector<std::complex> newrow(column);
        for (int j=0;j<column;j++)
            newrow[columnindex++]=data[i*column+j];
        back[rowindex++]=newrow;
    }
    return back;
}
//析构函数，不能使用析构函数的原因在于很多函数局部类型返回需要保留data数据，不能使局部类型析构造成接收对象不可用。
// 但是如果重写了赋值运算符，可以加上，因为在函数内，先赋值，在析构局部变量
CMatrix::~CMatrix(){
    row=0;column=0;
    if(data !=nullptr){
        delete[] data;
        data=nullptr;
    }
}