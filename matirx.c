#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include <complex.h>//用于复数的引用
#include "memory.h"

typedef unsigned int uint16_t;

struct complex
{
	float real;
	float image;
};

//struct complex complex;
typedef struct
{
    uint16_t row;
    uint16_t column;
    float **data;
	struct complex **complex_data;
}Matrix_t;



void clear_mat(Matrix_t* mat)
{
    uint16_t i, j;
    for (i = 0; i < mat->row; i++)
    {
        for (j = 0; j < mat->column; j++)
        {
            mat->data[i][j] = 0;
        }
    }
}
Matrix_t create_mat(uint16_t row, uint16_t column)
{
    Matrix_t mat;
    if (row <= 0||column<=0)
    {
        printf("error, in create_mat: row <= 0||column<=0\n");
        exit(1);
    }
    if (row > 0 && column > 0)
    {
        mat.row = row;
        mat.column = column;
        mat.data = (float **)malloc(row*sizeof(float *));//先指针的指针
		mat.complex_data=(complex **)malloc(row*sizeof(complex *));//先指针的指针
        if (mat.data == NULL)
        {
            printf("error, in create_mat: mat.data==NULL");
            exit(1);
        }
        uint16_t i;
        for (i = 0; i < row; i++)
        {
            *(mat.data + i) = (float *)malloc(column*sizeof(float));//再分配每行的指针
			*(mat.complex_data + i) = (complex *)malloc(column*sizeof(complex));//再分配每行的指针
            if (mat.data[i] == NULL)
            {
              printf("error, in create_mat: mat.data==NULL");
              exit(1);
            }
        }
        clear_mat(&mat);
    }
    return mat;
}


void free_mat(Matrix_t *mat)
{
    uint16_t i;
    for (i = 0; i < mat->row; i++)
	{
		free(mat->data[i]);/*释放行*/
		free(mat->complex_data[i]);/*释放行*/
	}
        
	    
    free(mat->data);/*释放头指针*/
	free(mat->complex_data);/*释放头指针*/
}

Matrix_t set_mat_data(Matrix_t mat,const float *data)
{
    uint16_t i, j;
    for (i = 0; i < mat.row; i++)
    {
        for (j = 0; j < mat.column; j++)
        {
            mat.data[i][j] = data[i*mat.column+j];
			printf("%f ",mat.data[i][j]);
        }
		printf("\n");
    }
	return mat;
}



complex add_complex(complex c1,complex c2)
{
	complex c;
	c.real=c1.real+c2.real;
	c.image=c1.image+c2.image;
	return c;
}

complex sub_complex(complex c1,complex c2)
{
	complex c;
	c.real=c1.real-c2.real;
	c.image=c1.image-c2.image;
	return c;
}

complex mut_complex(complex c1,complex c2)
{
	complex c;
	c.real=c1.real*c2.real-c1.image*c2.image;
    c.image=c1.real*c2.image+c1.image*c2.real;
	return c;
}

double abs_complex(complex c1)
{
	double c ;
	c=sqrt(c1.real*c1.real+c1.image*c1.image);
	return c;
}
complex conj_complex(complex c1)
{
	complex c;
	c.real=c.real;
	c.image=(-c1.image);
	return c;
}

Matrix_t mut_complex_matrix(Matrix_t mat1,Matrix_t mat2)
{
	
	Matrix_t mat=create_mat(mat1.row,mat1.column);
	for(int i=0;i<mat1.row;i++)
	{
		for(int j=0;j<mat1.column;j++)
		{
			mat.complex_data[i][j]=mut_complex(mat1.complex_data[i][j],mat2.complex_data[i][j]);
		}
	}
	return mat;
}
Matrix_t abs_complex_matrix(Matrix_t mat1)
{
	
	Matrix_t mat=create_mat(mat1.row,mat1.column);
	for(int i=0;i<mat1.row;i++)
	{
		for(int j=0;j<mat1.column;j++)
		{
			mat.data[i][j]=abs_complex(mat1.complex_data[i][j]);//复数的绝对值之后变为实数
		}
	}
	return mat;
}

Matrix_t real(Matrix_t mat)
{
	for(int i=0;i<mat.row;i++)
	{
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=mat.complex_data[i][j].real;
			
			//printf("%f+%fi\n", creal(complex.data[i][j]),cimag(complex.data[i][j]));
		}
		 printf("\n");
	}
	return mat;
}
Matrix_t imag(Matrix_t mat)
{
	for(int i=0;i<mat.row;i++)
	{
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=mat.complex_data[i][j].image;
			
			//printf("%f+%fi\n", creal(complex.data[i][j]),cimag(complex.data[i][j]));
		}
		 printf("\n");
	}
	return mat;
}
Matrix_t abs_matrix(Matrix_t mat1)//实数的绝对值
{
	
	Matrix_t mat=create_mat(mat1.row,mat1.column);
	for(int i=0;i<mat1.row;i++)
	{
		for(int j=0;j<mat1.column;j++)
		{
			//mat.data[i][j]=sqrt(mat1.complex_data[i][j].real*mat1.complex_data[i][j].real+mat1.complex_data[i][j].image*mat1.complex_data[i][j].image);//abs(mat1.complex_data[i][j]);
	     	mat.data[i][j]=sqrt(mat.data[i][j]);
		    
		}
	}
	return mat;
}

Matrix_t conj_matrix(Matrix_t mat1)
{
	Matrix_t mat=create_mat(mat1.row,mat1.column);
	for(int i=0;i<mat1.row;i++)
	{
		for(int j=0;j<mat1.column;j++)
		{
			mat.complex_data[i][j]=conj_complex(mat1.complex_data[i][j]);
		}
	}
	return mat;
	
}

complex trace_complex(Matrix_t mat)
{
	complex sum;
	for (int i = 0; i < mat.row; i++)
	{
		for(int j=0;j<mat.column;j++)
		{
			if(i==j)
			{
				sum=add_complex(sum,mat.complex_data[i][j]);
			}
		}

	}
	return sum;
	
}
Matrix_t randi(int iMin,int iMax,uint16_t m,uint16_t n)
{
     //initial(T,m,n); 
	 Matrix_t mat=create_mat(m,n);
     srand((unsigned int) time(NULL));
	 for(int i=0;i<mat.row;i++)
	 {
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=rand()%(iMax-iMin+1)+iMin;
			
			printf("%f ",mat.data[i][j]);
		}
		 printf("\n");
	 }
	return mat;
}

Matrix_t matrix_divide_num(Matrix_t mat,int num)//矩阵元素除以某个数
{
	
	for(int i=0;i<mat.row;i++)
	 {
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=mat.data[i][j]/num;
			
			printf("%f ",mat.data[i][j]);
		}
		 printf("\n");
	 }
	return mat;
}

Matrix_t matrix_multiply_num(Matrix_t mat,int num)//矩阵元素乘以某个数
{
	
	for(int i=0;i<mat.row;i++)
	 {
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=mat.data[i][j]*num;
			
			printf("%f ",mat.data[i][j]);
		}
		 printf("\n");
	 }
	return mat;
}

Matrix_t matrix_abs(Matrix_t mat)//矩阵元素乘以某个数
{
	Matrix_t mat_abs=create_mat(mat.row,mat.column);;
	for(int i=0;i<mat.row;i++)
	 {
	 	for(int j=0;j<mat.column;j++)
		{
			mat_abs.data[i][j]=abs(mat.data[i][j]);
			
			printf("%f ",mat_abs.data[i][j]);
		}
		 printf("\n");
	 }
	return mat_abs;
}
Matrix_t matrix_power_num(Matrix_t mat,int num)//矩阵元素求幂次方
{
	
	for(int i=0;i<mat.row;i++)
	 {
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=pow(mat.data[i][j],num);
			
			printf("%f ",mat.data[i][j]);
		}
		 printf("\n");
	 }
	return mat;
}
Matrix_t matrix_sqrt_num(Matrix_t mat)//开方
{
	
	for(int i=0;i<mat.row;i++)
	 {
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=sqrt(mat.data[i][j]);
			
			printf("%f ",mat.data[i][j]);
		}
		 printf("\n");
	 }
	return mat;
}
Matrix_t matrix_mean(Matrix_t mat)//求平均
{
	//free_mat(&mat);
	Matrix_t mean_mat;
	float sum;
	printf("mean:");
	for(int i=0;i<mat.column;i++)
	 {
	 	//mean_mat.data[0][i]=0;
		float sum = 0;
		for(int j=0;j<mat.row;j++)
		{
			sum+=mat.data[j][i];
		}
		mean_mat.data[0][i]=sum/mat.row;
		printf("%f ",mean_mat.data[0][i]);
	 }
	return mean_mat;
}

double matrix_mean_1(Matrix_t mat)//求平均,列数为1
{
	//free_mat(&mat);
	//Matrix_t mean_mat;
	float sum;
	printf("mean:");
	for(int i=0;i<1;i++)
	 {
	 	//mean_mat.data[0][i]=0;
		float sum = 0;
		for(int j=0;j<mat.row;j++)
		{
			sum+=mat.data[j][i];
		}
		sum=sum/mat.row;
		printf("%f ",sum);
	 }
	return sum;
}

#define PI 3.141592654
double gaussrand( )//获取正态分布随机数
{
    static double U, V;
    static int phase = 0;
    double Z;

    if(phase == 0)
    {
         U = rand() / (RAND_MAX + 1.0);
         V = rand() / (RAND_MAX + 1.0);
         Z = sqrt(-2.0 * log(U))* sin(2.0 * PI * V);
    }
    else
    {
         Z = sqrt(-2.0 * log(U)) * cos(2.0 * PI * V);
    }

    phase = 1 - phase;
    return Z;
}

Matrix_t randn(uint16_t m,uint16_t n)//m*n随机正态分布
{
     //initial(T,m,n); 
	 Matrix_t mat=create_mat(m,n);
     srand((unsigned int) time(NULL));
	 for(int i=0;i<mat.row;i++)
	 {
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=gaussrand();
			
			printf("%f ",mat.data[i][j]);
		}
		 printf("\n");
	 }
	return mat;
}

Matrix_t rand0_1(uint16_t m,uint16_t n)//0到1之间均匀伪随机数矩阵
{
	 Matrix_t mat=create_mat(m,n);
     srand((unsigned int) time(NULL));
	 for(int i=0;i<mat.row;i++)
	 {
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=(double)rand()/RAND_MAX;
			printf("%f ",mat.data[i][j]);
		}
		 printf("\n");
	 }
	return mat;
}




Matrix_t ones(uint16_t m,uint16_t n)//1矩阵
{
	Matrix_t mat=create_mat(m,n);
	for(int i=0;i<mat.row;i++)
	 {
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=1;
			
			printf("%d ",(unsigned char)(mat.data[i][j]));
		}
		 printf("\n");
	 }
	return mat;
}

Matrix_t cumsum(Matrix_t mat_t)//列向量求和
{
	Matrix_t mat;
	mat=mat_t;
	Matrix_t cumsum=create_mat(mat.row,mat.column+1);
	printf("row:%d column:%d\n",cumsum.row,cumsum.column);
	for(int i=0;i<cumsum.row;i++)
	 {
	 	for(int j=0;j<cumsum.column;j++)
		{
			
		    if(j==0)
			{
				cumsum.data[i][j]=0;
			}else if(j==1)
			{
				cumsum.data[i][j]=mat.data[i][0];
			}else
			{
				for(int k=0;k<j;k++)
				{
					cumsum.data[i][j]+=mat.data[i][k];
				}
			}
			
			
			printf("%d ",(unsigned char)(cumsum.data[i][j]));
		}
		 printf("\n");
	 }
	return cumsum;
}

int comp_ascend(const void *p1,const void *p2)//由小到大排序
{
    return  *((int*)p2)>*((int*)p1)?-1:1;
}
int comp_deascend(const void *p1,const void *p2)//由大到小排序
{
    return  *((int*)p2)>*((int*)p1)?1:-1;
}
Matrix_t sort(Matrix_t mat,uint16_t flag)//排序
{
	if(flag==1)//由小到大排序
	{
	   	for(int i=0;i<mat.row;i++)
	   {

			qsort(mat.data[i], sizeof(mat.column)+1,sizeof(float),comp_ascend);//调用函数qsort
			
			for(int j=0;j<mat.column;j++)
			{
				printf("%d ",(unsigned char)(mat.data[i][j]));
			}
	   }
	}else//由大到小排序
	{
		for(int i=0;i<mat.row;i++)
	   {

			qsort(mat.data[i], sizeof(mat.column),sizeof(float),comp_deascend);//调用函数qsort
			for(int j=0;j<mat.column;j++)
			{
				printf("%d ",(unsigned char)(mat.data[i][j]));
			}
	 }
	}
	return mat;
}

Matrix_t repmat(Matrix_t mat,uint16_t m,uint16_t n)//复制矩阵
{
	Matrix_t repmat=create_mat(mat.row*m,mat.column*n);
	for(int q=0;q<m;q++)
	{
		for(int i=0;i<mat.row;i++)
		 {
			for(int j=0;j<n;j++)
			{
				for(int k=0;k<mat.column;k++)
				{
					repmat.data[q][j]=mat.data[i][k];
					printf("%d ",(unsigned char)(repmat.data[q][j]));
				}			
			}
			printf("\n");
		 }
	}
	return repmat;
}

Matrix_t transpose_mat(Matrix_t mat)//矩阵转置
{
    Matrix_t mat_T;
    mat_T = create_mat(mat.column, mat.row);
    for (int i = 0; i < mat.column; i++)
    {
        for (int j = 0; j < mat.row; j++)
        {
            mat_T.data[i][j] = mat.data[j][i];
			printf("%d ",(unsigned char)(mat_T.data[i][j]));
        }
		
		printf("\n");
    }
    return mat_T;
}

Matrix_t sum(Matrix_t mat,uint16_t dim)//矩阵按行向量或列向量相加
{
    Matrix_t mat_S;
	if(dim == 1)//行向量元素相加
	{
		mat_S = create_mat(1,mat.column);

		for(int i=0;i<mat.column;i++)
		{
			for (int j = 0; j < mat.row; j++)
			{
				mat_S.data[0][i] += mat.data[j][i];
			}
			printf("%d ",(unsigned char)(mat_S.data[0][i]));
		}
		
		printf("\n");
	}
	else    //默认列向量相加 dim==2
	{
		mat_S = create_mat(mat.row,1);
		for (int i = 0; i < mat.row; i++)
		{
			for(int j=0;j<mat.column;j++)
			{
				mat_S.data[i][0] += mat.data[i][j];
			}
			printf("%d ",(unsigned char)(mat_S.data[i][0]));
			printf("\n");
		}
	}
    return mat_S;
}
Matrix_t mult_mat(Matrix_t mat1, Matrix_t mat2)
{
    Matrix_t mat;
    if (mat1.column != mat2.row)
    {
        printf("error,In mult_mat: mat1.column != mat2.row\n");
        exit(1);
    }   
    else
    {
        mat = create_mat(mat1.row, mat2.column);
        clear_mat(&mat);
        for (int i = 0; i < mat1.row; i++)
        {
            for (int j = 0; j < mat2.column; j++)
            {
                for (int m = 0; m < mat1.column; m++)
                {
                    mat.data[i][j] += mat1.data[i][m] * mat2.data[m][j];
                }
				printf("%d ",(unsigned char)(mat.data[i][j]));
            }
			printf("\n");
        }
    }
    return mat;
}
Matrix_t divide_mat(Matrix_t mat1, Matrix_t mat2)
{
    Matrix_t mat;
    if (mat1.column != mat2.row)
    {
        printf("error,In mult_mat: mat1.column != mat2.row\n");
        exit(1);
    }   
    else
    {
        mat = create_mat(mat1.row, mat2.column);
        clear_mat(&mat);
        for (int i = 0; i < mat1.row; i++)
        {
            for (int j = 0; j < mat2.column; j++)
            {
                for (int m = 0; m < mat1.column; m++)
                {
                    mat.data[i][j] += mat1.data[i][m] / mat2.data[m][j];
                }
				printf("%d ",(unsigned char)(mat.data[i][j]));
            }
			printf("\n");
        }
    }
    return mat;
}

Matrix_t add_sub_mat(Matrix_t mat1, Matrix_t mat2,int k)//矩阵相加减
{
    Matrix_t mat;
   /* if (mat1.column != mat2.row)
    {
        printf("error,In mult_mat: mat1.column != mat2.row\n");
        exit(1);
    }   
    else*/
    {
        mat = create_mat(mat1.row, mat1.column);

        for (int i = 0; i < mat1.row; i++)
        {
            for (int j = 0; j < mat1.column; j++)
            {
				mat.data[i][j] = mat1.data[i][j] +((k>0)?mat2.data[i][j]:-mat2.data[i][j]);
            }
			
        }
    }
    return mat;
}

unsigned length(Matrix_t mat)
{
	return mat.row>=mat.column?mat.row:mat.column;
}
unsigned size(Matrix_t mat,unsigned mode)//返回矩阵行数或列数
{
	if(mode==1)
	{
		return mat.row;
	}else
	{
		return mat.column;
	}
}

double trace(Matrix_t mat)
{
	double sum;
	for (int i = 0; i < mat.row; i++)
	{
		for(int j=0;j<mat.column;j++)
		{
			if(i==j)
			{
				sum+=mat.data[i][j];
			}
		}

	}
	return sum;
	
}
Matrix_t mat_spdiags(Matrix_t mat,int k,int m,int n)
{
	Matrix_t mat_spdiags;
    mat_spdiags=create_mat(m,n);
	for(int i=0;i<m;i++)
    {
		for(int j=0;j<n;j++)
		{
			if(i==j)
			{
				mat_spdiags.data[i][j]=mat.data[j][0];
			}
			else
			{
				mat_spdiags.data[i][j]=0;
			}
			printf("%f ",mat_spdiags.data[i][j]);
		}
		printf("\n");
    }
	return mat_spdiags;
        
}

/*Matrix_t real(Matrix_t mat)
{
	for(int i=0;i<mat.row;i++)
	{
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=creal(mat.data[i][j]);
			
			//printf("%f+%fi\n", creal(complex.data[i][j]),cimag(complex.data[i][j]));
		}
		 printf("\n");
	}
	return mat;
}
Matrix_t imag(Matrix_t mat)
{
	for(int i=0;i<mat.row;i++)
	{
	 	for(int j=0;j<mat.column;j++)
		{
			mat.data[i][j]=cimag(mat.data[i][j]);
			
			//printf("%f+%fi\n", creal(complex.data[i][j]),cimag(complex.data[i][j]));
		}
		 printf("\n");
	}
	return mat;
}*/
/*Matrix_t conj(Matrix_t mat_real,Matrix_t mat_imag)
{
	Matrix_t mat_conj;
    mat_conj=create_mat(mat_real.row,mat_real.column);
	for(int i=0;i<mat.row;i++)
	{
	 	for(int j=0;j<mat.column;j++)
		{
			mat_conj.data[i][j]=(mat_real.data[i][j]- mat_imag.data[i][j]);
			
			//printf("%f+%fi\n", creal(complex.data[i][j]),cimag(complex.data[i][j]));
		}
		 printf("\n");
	}
	return mat_conj;
}*/

double norm_sub(Matrix_t mat1,Matrix_t mat2)//两个行向量的模相减取模，范数
{
	double sum;
	for (int i = 0; i < mat1.row; i++)
	{
		for(int j=0;j<mat1.column;j++)
		{
			sum += pow(mat2.data[i][j]-mat1.data[i][j],2);
		}
	}
	return sqrt(sum);
}

double norm(Matrix_t mat)//行向量的取模，范数
{
	double sum;
	for (int i = 0; i < mat.row; i++)
	{
		for(int j=0;j<mat.column;j++)
		{
			sum += pow(mat.data[i][j],2);
		}
	}
	
	return sqrt(sum);
}
Matrix_t findpeaks(Matrix_t mat,int row)//按照matlab得出结果后要sortstr排序，下降;row=0:全部行数据返回，1:返回第一行，2返回第二行
{
	Matrix_t findpeaks;
	Matrix_t findpeaks1;
	int n;
	findpeaks = create_mat(2, mat.column);//第一行为pks，第二行为locs
	for (int i = 0; i < mat.row; i++)
	{
		for(int j=0;j<(mat.column-2);j++)
		{
			if((mat.data[i][j]<=mat.data[i][j+1]) && (mat.data[i][j+1]>=mat.data[i][j+2]))
			{
				if(mat.data[i][j+1] != mat.data[i][j])
				{
					findpeaks.data[i][n]=mat.data[i][j+1];
					findpeaks.data[i+1][n]=j+2;
					printf("%f  %f",findpeaks.data[i][n],findpeaks.data[i+1][n]);
					n++;
				}

			}
		}
	}
	findpeaks1 = create_mat(2, n-1);//第一行为pks，第二行为locs
	printf("\n");
	if(row==0)
	{
		for (int i = 0; i < 2; i++)
		{
			for(int j=0;j<n-1;j++)
			{
				findpeaks1.data[i][j]=findpeaks.data[i][j+1];
				printf("%f  ",findpeaks1.data[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		//findpeaks1=add_sub_mat(findpeaks,2,k);
		//
		return sort(findpeaks1,2);//从大到小排序
	}else if(row==1)
	{
		
		for (int i = 0; i < 1; i++)
		{
			for(int j=0;j<n-1;j++)
			{
				findpeaks1.data[i][j]=findpeaks.data[i][j+1];
				printf("%f  ",findpeaks1.data[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		//findpeaks1=add_sub_mat(findpeaks,2,k);
		//
		return sort(findpeaks1,2);//从大到小排序
	}else if(row==2)
	{
		for (int i = 1; i < 2; i++)
		{
			for(int j=0;j<n-1;j++)
			{
				findpeaks1.data[i][j]=findpeaks.data[i][j+1];
				printf("%f  ",findpeaks1.data[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		//findpeaks1=add_sub_mat(findpeaks,2,k);
		//
		return sort(findpeaks1,2);//从大到小排序
	}
	
}

Matrix_t eye(unsigned num)//对角线全为1，其余为0
{
	Matrix_t eye;
	eye = create_mat(num, num);
	for (int i = 0; i < num; i++)
	{
		for(int j=0;j<num;j++)
		{
			if(i==j)
			{
				eye.data[i][j]=1;
			}else
			{
				eye.data[i][j]=0;
			}
		}
	}
	return eye;
}

void swap(Matrix_t input,Matrix_t output,int dimension,int first_row,int second_row)
{
	double temp_row1[dimension];
	double temp_row2[dimension];
    int i;
    for(i=0;i<dimension;i++)
    {
        temp_row1[i]=input.data[first_row][i];
        temp_row2[i]=output.data[first_row][i];
    }
    for(i=0;i<dimension;i++)
    {
        input.data[first_row][i]=input.data[second_row][i];
        output.data[first_row][i]=output.data[second_row][i];
        input.data[second_row][i]=temp_row1[i];
        output.data[second_row][i]=temp_row2[i];
    }
}

void reorderOutput(Matrix_t output,int dimension)
{
	Matrix_t temp;
	temp=create_mat(dimension,dimension);
	for(int i=1;i<dimension;++i)
	{
		memcpy(temp.data[i-1],output.data[i],sizeof(double)*dimension);
	}
	memcpy(temp.data[dimension-1],output.data[0],sizeof(double)*dimension);
	for(int i=0;i<dimension;++i)
        memcpy(output.data[i],temp.data[i],sizeof(double)*dimension);
}
//高斯消元法求矩阵逆
Matrix_t inverse_Matrix(Matrix_t input,int dimension)
{
    Matrix_t output=create_mat(dimension,dimension);
	int i,j,k;    
	output=eye(dimension);//将输出矩阵初始化为单位矩阵
    for(i=0;i<dimension;++i)  //依次处理每一列
    {
        for(j=0;i<dimension;++j)  //如果当前行当前列值为0，做行变换
        {
            if(input.data[j][i]!=0)
            {
                swap(input,output,dimension,0,j);
                break;
            }
        }
        for(j=0;j<dimension;++j)  //依次处理每一行
        {
            if(j==0)  //如果是第一行，将input[j][i]设置为1，其他元素均除以input[i][i]
            {
                for(k=dimension-1;k>=0;--k)
                    output.data[j][k]/=input.data[j][i];
                for(k=dimension-1;k>=i;--k)
                    input.data[j][k]/=input.data[j][i];
            }
            else  //如果不是第一行，将每一行的input[j][i]设置为0，该行其他元素都要倍数的减去第一行对应位置上的值
            {
                for(k=dimension-1;k>=0;--k)
                    output.data[j][k]=output.data[j][k]-input.data[j][i]/input.data[0][i]*output.data[0][k];
                for(k=dimension-1;k>=i;--k)
                    input.data[j][k]=input.data[j][k]-input.data[j][i]/input.data[0][i]*input.data[0][k];
            }
        }
        swap(input,output,dimension,0,(i+1)%dimension);  //每次都将下一次需要处理的行和当前的第一行交换
    }
	
    reorderOutput(output,dimension); //因为之前的交换操作，行顺序乱了，需要重新排列一下，即把第一行的数据放到最后一行后面
    return output;
}
/*Matrix_t complex(Matrix_t mat)
{
	
	Matrix_t complex;
	//complex double complex.data;
	complex=create_mat(mat.row,mat.column);
	//mat = randn(m,n);
	for(int i=0;i<mat.row;i++)
	{
	 	for(int j=0;j<mat.column;j++)
		{
			complex.data[i][j]=(mat.data[i][j]+ mat.data[i][j]);
			
			//printf("%f+%fi\n", creal(complex.data[i][j]),cimag(complex.data[i][j]));
		}
		 printf("\n");
	}
	return complex;

}*/

Matrix_t arrange_row(Matrix_t mat)//所有的元素按行排列,相当于matlab中的Y(:)函数
{
	Matrix_t arrange_row=create_mat(mat.row*mat.column,1);
	for(int i=0;i<mat.column;i++)
	{
		for(int j=0;j<mat.row;j++)
		{
			arrange_row.data[i*mat.row+j][0]=mat.data[j][i];
			printf("%f \n",arrange_row.data[i*mat.row+j][0]);
		}
	}
	return 	arrange_row;
}

Matrix_t index_value(Matrix_t mat,Matrix_t index,int start,int end)  //根据start-end索引找到相应的值
{
	printf("index_value\n");
	if(mat.row==1)
	{
		Matrix_t index_value;
		index_value=create_mat(start,end);
		for(int i=0;i<end;i++)
		{
			for(int j=0;j<mat.column;j++)
			{
				if(j==index.data[1][i])
				{
					index_value.data[0][i]=mat.data[0][j];
					printf("%f ",index_value.data[0][i]);
					break;
				}
			}
		}
		return index_value;

	}else if(mat.column==1)
	{
		Matrix_t index_value=create_mat(end,start);
		for(int i=0;i<end;i++)
		{
			for(int j=0;j<mat.row;j++)
			{
				if(j==index.data[1][i])
				{
					index_value.data[i][0]=mat.data[j][0];
					printf("%f \n",index_value.data[i][0]);
					break;
				}
			}
		}
		
		return index_value;
	}
}

Matrix_t index_change(Matrix_t mat1,Matrix_t mat2)//比对两个矩阵值，返回比对的索引
{
	Matrix_t index_change=create_mat(1,mat2.row*mat2.column);
	for(int i=0;i<mat2.row;i++)
	{
		for(int j=0;j<mat2.column;j++)
		{
			for(int k=0;k<mat1.row;k++)
	        {
		         for(int t=0;t<mat1.column;t++)
	             {
					if(mat2.data[i][j]==mat1.data[k][t])
					{
						index_change.data[0][i*mat2.column+j]=k*mat1.column+t;
						break;
					}
				 }
            }
		}
	}
	return index_change;
}
Matrix_t index_sort(Matrix_t mat1,Matrix_t mat2)//根据mat2的索引，排序mat1的顺序
{
	Matrix_t mat=create_mat(mat1.row,mat1.column);
	printf("\n重新排序：");
	for(int i=0;i<mat2.column;i++)
	{
		mat.data[0][i]=mat1.data[0][(int)mat2.data[0][i]];
		printf("%f ",mat.data[0][i]);
	}
	return mat;
}

void fun_SAM3Res(Matrix_t Y,Matrix_t A,Matrix_t DAS_init,Matrix_t DOAscan,Matrix_t DOA)//
	   //fun_SAM3Res(Y,A,DAS_init,DOAscan,DOA)
{
	double Numsources =length(DOA);
	double threshold=pow(10,-6);
	double maxIter=30;
	double M=size(A,1);//返回行数
	double thetaNum=size(A,2);//返回列数
	double t_samples=size(Y,2);//返回列数
	Matrix_t R_N;
	R_N = matrix_divide_num(mult_mat(Y,transpose_mat(Y)),t_samples);//原函数R_N = (Y*Y')/t_samples;
	//sigma = mean(abs(Y(:)).^2); 
	//先求Y(:),abs(),.^2,再求mean
	//matrix_mean(matrix_power_num())
	double sigma;
	sigma=matrix_mean_1(matrix_power_num(abs_matrix(arrange_row(Y)),2));//原函数sigma = mean(abs(Y(:)).^2);
	Matrix_t p_vec_Old;
	p_vec_Old =matrix_power_num(abs_matrix(DAS_init),2);//原函数p_vec_Old = abs(DAS_init).^2;
	Matrix_t R;
	Matrix_t Rinv;
	Matrix_t Rinv_A;
	Matrix_t diag_A_Rinv_A;
	Matrix_t mut_tmp,tmp;
	Matrix_t tmp_divide,p_vec;
	Matrix_t tmp_trace;
	double real1,real2;
	double norm1,norm2,p_diffs_ratio;
	
	for(int iterIdx=1;iterIdx<maxIter;iterIdx++)//原函数for iterIdx = 1:maxIter
	{
		//R=A*mat_spdiags(p_vec_Old,0,thetaNum,thetaNum)*transpose_mat(Y)+sigma*eye(M);//原函数 R =  A*spdiags(p_vec_Old, 0, thetaNum, thetaNum)*A' + sigma*eye(M);
		Matrix_t mat_spdiags_1,transpose_mat_1,eye_1,a1,b1,c1,d1;
        mat_spdiags_1=mat_spdiags(p_vec_Old,0,thetaNum,thetaNum);
        a1=mult_mat(A,mat_spdiags_1);
        transpose_mat_1=transpose_mat(Y);
        b1=mult_mat(a1,transpose_mat_1);
        eye_1=eye(M);
        c1=matrix_multiply_num(eye_1,sigma);
        d1=add_sub_mat(b1,c1,1);
        
		Rinv=inverse_Matrix(R,R.row);//原函数Rinv=inv(R)
		Rinv_A=	mult_mat(Rinv,A);//原函数Rinv_A = Rinv*A;
		diag_A_Rinv_A=transpose_mat(sum(mult_mat(conj_matrix(A),Rinv_A),1));//原函数diag_A_Rinv_A = sum(conj(A).*Rinv_A, 1).';
		mut_tmp=mult_mat(mult_mat(mult_mat(Rinv,R_N),Rinv),A);
		tmp=transpose_mat(sum(mult_mat(conj_matrix(A),mut_tmp),1));//原函数tmp = sum(conj(A).* (Rinv*R_N*Rinv*A), 1).';
		tmp_divide=divide_mat(tmp,diag_A_Rinv_A);
		p_vec=mult_mat(p_vec_Old,tmp_divide);//原函数p_vec = p_vec_Old.*(tmp./diag_A_Rinv_A);
		tmp_trace=mult_mat(mult_mat(Rinv,Rinv),R_N);
		
		complex trace1,trace2;
		trace1=trace_complex(tmp_trace);
		trace2=trace_complex(mult_mat(Rinv,Rinv));
		real1=trace1.real;//real(trace1);
		real2=trace2.real;//real(trace2);
		sigma=real1/real2;//原函数sigma = real(trace(Rinv*Rinv*R_N))/real(trace(Rinv*Rinv));
		
		norm1=norm_sub(p_vec_Old,p_vec);
		norm2=norm(p_vec_Old);
		p_diffs_ratio=norm1/norm2;//原函数p_diffs_ratio = norm(p_vec_Old-p_vec)/norm(p_vec_Old);
		if(p_diffs_ratio<threshold)//原函数if p_diffs_ratio < threshold
			break;
		p_vec_Old = p_vec;
	}
	
	p_vec = real(p_vec);
	
	Matrix_t tmp_findpeaks,pks,index;
	//tmp_findpeaks=findpeaks(p_vec);
	pks=findpeaks(p_vec,1);
	index=findpeaks(p_vec,2);//原函数[pks index]=findpeaks(p_vec, 'sortstr', 'descend');
	
	double normal,noisepower;
	Matrix_t Distance,Detected_powers;
	if(length(index)<Numsources)
	{
    	normal = 0;
    	//Distance = 0;//NaN
     	//p_vec = 0;
    	//Detected_powers = 0;
    	noisepower = sigma;
	}
	Matrix_t Detected_DOAs;
	Detected_DOAs=index_value(DOAscan,index,1,Numsources);//原函数Detected_DOAs = DOAscan(index(1:Numsources));
	Matrix_t Detected_DOAs_tmp;
	Detected_DOAs_tmp=sort(Detected_DOAs,1);//1为按照上升排序： 原函数 
	Matrix_t IXsort;
	IXsort=index_change(Detected_DOAs,Detected_DOAs_tmp);   //[Detected_DOAs, IXsort] = sort(Detected_DOAs, 'ascend');
	Detected_DOAs=Detected_DOAs_tmp;
	//Matrix_t Distance;
	Distance=add_sub_mat(Detected_DOAs, DOA,-1);//矩阵相加减//Distance = Detected_DOAs - DOA;
	normal = 1;
	//Matrix_t Detected_powers;
	Detected_powers=index_value(p_vec,index,1,Numsources);//Detected_powers = p_vec(index(1:Numsources));
	Detected_powers=index_sort(Detected_powers,IXsort);//Detected_powers =  Detected_powers(IXsort);//还不知到怎么获取索引
	noisepower = sigma;//
}

int main(void)
{

	//complex double noisenew;
	Matrix_t mat2;
	Matrix_t mat;
    //answer = exp (10);
    //printf("e^10 =%f\n", answer);
	//mat=randi(-5,5,6,3);
	//mat1=randi(1,5,2,3);
	mat2=randi(1,5,4,5);
	//mat=randn(3,3);

    //mat=rand0_1(3,4);
	//mat=sum(mat,2);
	//mat=mult_mat(mat1,mat2);
	//mat=matrix_multiply_num(mat,2);
	//mat=matrix_abs(mat);
	//answer=trace(mat);
	//printf("answer=%f",answer);
	//mat=mat_spdiags(mat2,0,3,3);
	return 0;
}