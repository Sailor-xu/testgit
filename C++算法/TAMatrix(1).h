//
//  TAMatrix.hpp
//  
//
//  Created by Hongren Chen on 27/09/2017.
//

#ifndef TAMATRIX_H
#define TAMATRIX_H

#include <iostream>
#include <string.h>
#include <cmath>
#include <complex>
#include <vector>

class Matrix {
    int rows;
    int columns;
    int size;
    double *p;
public:
    Matrix();
    Matrix(int, int , double*);
    Matrix(int, int);
    Matrix(int, int, double);
    Matrix(const Matrix&);
    ~Matrix();
    
    const Matrix& operator+() const;
    const Matrix operator-() const;
    Matrix& operator++();
    const Matrix operator++(int);
    Matrix& operator--();
    const Matrix operator--(int);
    
    const Matrix operator+(const Matrix&) const;
    const Matrix operator-(const Matrix&) const;
    const Matrix operator*(const Matrix&) const;
    
    const Matrix operator+(const double&) const;
    const Matrix operator-(const double&) const;
    const Matrix operator*(const double&) const;
    const Matrix operator/(const double&) const;
    
    friend const Matrix operator+(const double&, const Matrix&);
    friend const Matrix operator-(const double&, const Matrix&);
    friend const Matrix operator*(const double&, const Matrix&);
    
    Matrix& operator=(const Matrix&);
    Matrix& operator+=(const Matrix&);
    Matrix& operator-=(const Matrix&);
    Matrix& operator*=(const Matrix&);
    
    Matrix& operator+=(const double&);
    Matrix& operator-=(const double&);
    Matrix& operator*=(const double&);
    Matrix& operator/=(const double&);
    
    bool operator==(const Matrix&) const;
    friend std::ostream& operator<<(std::ostream&, const Matrix&);
    
    double& takeElem(int, int) const;
    Matrix takeRow(int) const;
    Matrix takeColumn(int) const;
    Matrix takeSubMatrix(int, int, int, int) const;

    int& getRows() const;
    int& getColumns() const;
    
    Matrix inverse() const;
    Matrix transpose() const;
    double detMatrix() const;


    //double norm2();
    //double mean();
    Matrix eye();
    Matrix diag();
    Matrix cov(bool flag = true);

    void QR(Matrix &, Matrix &)const;
    Matrix eig_val(unsigned _iters = 1000);
    Matrix eig_vect(unsigned _iters = 1000);
};

const double DOUBLEERROR = 1.0e-15;  //最小double，比较误差
//比较两double浮点数不相等
bool FloatNotEqual(double lhs, double rhs)
{
    if (std::abs(lhs - rhs) >= DOUBLEERROR)
        return true;
    else
        return false;
}

//定义复数矩阵类
class CMatrix{
//公共数据
public:
    std::complex* data;
    int row;
    int column;
    //构造函数，构建空矩阵
    CMatrix();
    //由Complex转化为1x1矩阵
    CMatrix(std::complex);
    //由实矩阵转化为复数矩阵
    CMatrix(Matrix& );
    //由实部和虚部组成复数矩阵
    CMatrix(Matrix&,Matrix&);

    //构建矩阵，不赋值
    CMatrix(int ,int );
    //构建值矩阵，矩阵中每个元素都等于给定值
    CMatrix(int ,int ,std::complex );
    //根据数组构建矩阵,复制元素
    CMatrix(int ,int ,double*,double*);
    //根据数组构建矩阵,复制元素
    CMatrix(double*,double*,int ,int );
    //根据数组构建矩阵,复制元素
    CMatrix(int ,int ,std::complex *);
    //根据数组构建矩阵,复制元素
    CMatrix(std::complex *,int ,int );
    //根据向量构建矩阵
    CMatrix(std::vector<std::complex> );
    //根据向量构建矩阵
    CMatrix(std::vector<std::complex> ,int ,int );
    //根据向量构建矩阵
    CMatrix(std::vector<std::vector<std::complex>> );
    //复制构造函数
    CMatrix(const CMatrix& );
    //重载赋值运算符，保留在matlab中赋值等于深复制的效果
    CMatrix& operator=(const CMatrix& );
    //重载()运算符,下标从0开始的习惯，读取时调用
    std::complex operator()(int ,int ) const;
    //重载()运算符,下标从0开始的习惯，写入时调用
    std::complex& operator()(int ,int );
    //重载()运算符,下标从0开始的习惯，读取时调用
    std::complex operator()(int ) const;
    //重载()运算符,下标从0开始的习惯
    std::complex& operator()(int );

    //重载加法运算符
    CMatrix operator+(std::complex );
    //重载加法运算符
    CMatrix& operator+=(std::complex );
    //重载加法运算符
    CMatrix operator+(CMatrix& );
    //重载加法运算符
    CMatrix& operator+=(CMatrix& );
    //重载减法运算符
    CMatrix operator-(std::complex );
    //重载减法运算符
    CMatrix& operator-=(std::complex );
    //重载减法运算符
    CMatrix operator-(CMatrix& );
    //重载减法运算符
    CMatrix& operator-=(CMatrix& );
    //重载乘法运算符
    CMatrix operator*(std::complex );
    //重载乘法运算符
    CMatrix& operator*=(std::complex );
    //重载乘法运算符
    CMatrix operator*(CMatrix& );
    //重载除法运算符
    CMatrix operator/(std::complex );
    //重载除法运算符
    CMatrix& operator/=(std::complex );
    //重载矩阵等号==运算符
    bool operator==(CMatrix& );
    //重载矩阵等号!=运算符
    bool operator!=(CMatrix& );
    //重载一元操作符：正号
    CMatrix operator +() const;
    //重载一元操作符：负号
    CMatrix operator -() const;
    //矩阵转置，不改变源矩阵
    CMatrix transpose();

    //矩阵左扩展，改变了源矩阵
    void append_left(CMatrix& );
    //向量首部添加元素
    void append_left(std::complex );
    //矩阵右扩展，改变了源矩阵
    void append_right(CMatrix& );
    //向量尾部添加元素
    void append_right(std::complex );
    //矩阵上扩展，改变了源矩阵
    void append_top(CMatrix& );
    //向量首部添加元素
    void append_top(std::complex );
    //矩阵下扩展，改变了源矩阵
    void append_bottom(CMatrix& );
    //向量尾部添加元素
    void append_bottom(std::complex );
    //删除指定行，改变了源矩阵
    void remove_row(int );  //删除指定行
    //删除指定行，改变了源矩阵
    void remove_row(int ,int );  //删除多行,包含两个下标
    //删除指定列，改变了源矩阵
    void remove_column(int );  //删除指定列
    //删除指定列，改变了源矩阵
    void remove_column(int ,int );  //删除多列,包含两个下标
    //删除指定的行和列
    void remove_row_column(int ,int );
    //删除多行和多列
    void remove_row_column(int ,int ,int ,int );
    //替换局部矩阵，改变了源矩阵，
    void replace(CMatrix &,int ,int );

    //点运算
    CMatrix dot(CMatrix& ,std::string );

    //查询是否是向量
    bool isVector();
    //实部矩阵
    Matrix real();
    //虚部矩阵
    Matrix imag();
    //模值矩阵
    Matrix abs();
    //角度矩阵
    Matrix angle();
    //共轭矩阵
    CMatrix conj();

    //按行优先转化为series
    std::complex* toSeries();
    //转化为二位数组
    std::complex** toArray();
    //按行优先转化为一维vector
    std::vector<std::complex> toVector();
    //按行优先转化为二维vector
    std::vector<std::vector<std::complex>> toVector2();
    //析构函数，不能使用析构函数的原因在于很多函数局部类型返回需要保留data数据，不能使局部类型析构造成接收对象不可用。
    // 但是如果重写了赋值运算符，可以加上，因为在函数内，先赋值，在析构局部变量
    ~CMatrix();

};

#endif /* TAMATRIX_H */
