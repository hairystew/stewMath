#pragma once

#define PI 3.1415926535897932384626433832795
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>


#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif





namespace sm {




    //-----------------------------GENERIC 3D VECTOR-----------------------------------------
    template <class T>
    struct vec3d_generic //Creating 3-D vector structure
    {
        T x = 0.0;
        T y = 0.0;
        T z = 0.0;

        vec3d_generic() : x(0), y(0), z(0) {}
        vec3d_generic(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}
        vec3d_generic(const vec3d_generic& v) : x(v.x), y(v.y), z(v.z) {}


        //Adding utility functions

        T mag() {
            return sqrt(x * x + y * y + z * z);
        }

        vec3d_generic norm() {
            T r = 1 / mag();
            return vec3d_generic(x * r, y * r, z * r);
        }


        //Using operator overloading to allow the vec3d structor to behave like a primative

        vec3d_generic operator + (const vec3d_generic& rightSide)
        {
            return vec3d_generic(this->x + rightSide.x, this->y + rightSide.y, this->z + rightSide.z);
        }

        vec3d_generic operator - (const vec3d_generic& rightSide)
        {
            return vec3d_generic(this->x - rightSide.x, this->y - rightSide.y, this->z - rightSide.z);
        }

        vec3d_generic& operator += (const vec3d_generic& rightSide)
        {
            this->x += rightSide.x; this->y += rightSide.y; this->z += rightSide.z;
            return *this;
        }

        vec3d_generic& operator -= (const vec3d_generic& rightSide)
        {
            this->x -= rightSide.x; this->y -= rightSide.y; this->z -= rightSide.z;
            return *this;
        }

        vec3d_generic& operator *= (const T rightSide) //Multiplying by a scalar
        {
            this->x -= rightSide; this->y -= rightSide; this->z -= rightSide;
            return *this;
        }

        vec3d_generic& operator /= (const T rightSide) //Dividing by a scalar
        {
            this->x -= rightSide; this->y -= rightSide; this->z -= rightSide;
            return *this;
        }

        vec3d_generic operator * (const T rightSide) //Multiplying by a scalar
        {
            return vec3d_generic(this->x * rightSide, this->y * rightSide, this->z * rightSide);
        }



        vec3d_generic operator / (const T rightSide) //Dividing by a scalar
        {
            return vec3d_generic(this->x / rightSide, this->y / rightSide, this->z / rightSide);
        }


    };

    typedef vec3d_generic<double> vec3d;

    // overloading "<<" operator to allow for ease of printing to console, i.e "std::cout << genericVector"
    // I believe it is good practice to overload input operators to the iostream library as a nonmember function, however I could be wrong

    template <class T>
    std::ostream& operator<<(std::ostream& os, const vec3d_generic<T>& vec) {
        os << "{" << vec.x << ", " << vec.y << ", " << vec.z << "}";
        return os;
    }

    // Creating funcitons for cross product and dot product
    template <class T>
    vec3d_generic<T> cross(const vec3d_generic<T> vec1, const vec3d_generic<T> vec2) {
        return{ (vec1.y * vec2.z - vec2.y * vec1.z), -(vec1.x * vec2.z - vec2.x * vec1.z), (vec1.x * vec2.y - vec2.x * vec1.y) };
    }
    template <class T>
    double dot(const vec3d_generic<T> vec1, const vec3d_generic<T> vec2) {
        return(vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z);
    }

    template <class T>
    vec3d_generic<T> rot(const vec3d_generic<T> vec, const double angle, const char rotAxis) // angle in degrees, use single quotes when defining rotAxis... i.e. 'x' 
    {
        switch (rotAxis) {
        case 'x':
            return { vec.x, vec.y * cos(angle * PI / 180) - vec.z * sin(angle * PI / 180), vec.y * sin(angle * PI / 180) + vec.z * cos(angle * PI / 180) };
        case 'y':
            return { vec.x * cos(angle * PI / 180) + vec.z * sin(angle * PI / 180) ,vec.y ,-vec.x * sin(angle * PI / 180) + vec.z * cos(angle * PI / 180) };
        case 'z':
            return { vec.x * cos(angle * PI / 180) - vec.y * sin(angle * PI / 180), vec.x * sin(angle * PI / 180) + vec.y * cos(angle * PI / 180), vec.z };
        }
    }



    // MATRIX OPERATIONS

    template <class T>
    class Matrix
    {

    private:
        T **data; //Values in the Matrix
    public:



        int n; //Rows
        int m; //Cols


        T& operator()(const unsigned& n, const unsigned& m) {
            return this->data[n][m];
        }

        const T& operator()(const unsigned& n, const unsigned& m) const{
            return this->data[n][m];
        }

        //Initialized a N x M matrix to _init (N rows and M cols)
        Matrix<T>(int _n, int _m, const T& _init) : n(_n), m(_m) {
            data = new T* [n];
            for (unsigned i = 0; i < n; ++i) {
                data[i] = new T[m];
                for (unsigned j = 0; j < m; j++) {
                    data[i][j] = _init;
                }
            }

        }

        //Initialized a square matrix to identity
        Matrix<T>(int _n) : n(_n), m(_n) {
            data = new T * [n];
            for (unsigned i = 0; i < n; ++i) {
                data[i] = new T[m];
                for (unsigned j = 0; j < m; j++) {
                    if (i == j) { data[i][j] = 1; }
                }
            }

        }

        //copy constructor
        Matrix<T>(const Matrix<T>& rhs) {
            n = rhs.n;
            m = rhs.m;
            data = new T * [n];
            for (unsigned i = 0; i < n; ++i) {
                data[i] = new T[m];
                for (unsigned j = 0; j < m; j++) {
                    data[i][j] = rhs.data[i][j];
                }
            }
        }

        ~Matrix<T>() {
            for (unsigned i = 0; i < n; ++i) {
                delete[] data[i];
            }
            delete[] data;
        }

        Matrix<T> operator=(const Matrix<T>& rhs) {
            if (&rhs == this)
                return *this;

            for (unsigned i = 0; i < n; ++i) {
                delete[] data[i];
            }
            delete[] data;

            data = new T * [n];
            for (unsigned i = 0; i < n; ++i) {
                data[i] = new T[m];
                for (unsigned j = 0; j < m; j++) {
                    data[i][j] = rhs(i, j);
                }
            }
            n = rhs.n;
            m = rhs.m;
            return *this;
        }



        const void print() const{
            for (unsigned i = 0; i < n; i++) {
                for (unsigned j = 0; j < m; j++) {
                    std::cout << data[i][j] << ", ";
                }
                std::cout << std::endl;
            }
        }

        void print() {
            for (unsigned i = 0; i < n; i++) {
                for (unsigned j = 0; j < m; j++) {
                    std::cout << data[i][j] << ", ";
                }
                std::cout << std::endl;
            }
        }


        Matrix<T> transpose() {
            ASSERT(n == m, "MATRIX NOT SQUARE");
            Matrix<T> res(n, m, 0);
            for (unsigned i = 0; i < n; i++) {
                for (unsigned j = 0; j < m; j++) {
                    res(i, j) = data[j][i];
                }
            }
            return res;
        }

        //matrix addition
        Matrix<T> operator+(const Matrix<T>& rhs) {
            ASSERT(rhs.n == n || rhs.m == m, "MATRIX SIZE MISMATCH");
            Matrix<T> res(n ,m , 0);
            for (unsigned i = 0; i < n; i++) {
                for (unsigned j = 0; j < m; j++) {
                    res(i, j) = rhs(i, j) + this->data[i][j];
                }
            }
            return res;
        }

        //scalar addition
        Matrix<T> operator+(const T& rhs) {
            ASSERT(rhs.n == n || rhs.m == m, "MATRIX SIZE MISMATCH");
            Matrix<T> res(n, m, 0);
            for (unsigned i = 0; i < n; i++) {
                for (unsigned j = 0; j < m; j++) {
                    res(i, j) = rhs + this->data[i][j];
                }
            }
            return res;
        }

        //matrix subtraction
        Matrix<T> operator-(const Matrix<T>& rhs) {
            ASSERT(rhs.n == n || rhs.m == m, "MATRIX SIZE MISMATCH");
            Matrix<T> res(n, m, 0);
            for (unsigned i = 0; i < n; i++) {
                for (unsigned j = 0; j < m; j++) {
                    res(i, j) = this->data[i][j] - rhs(i, j);
                }
            }
            return res;
        }

        //scalar subtraction
        Matrix<T> operator-(const T& rhs) {
            ASSERT(rhs.n == n || rhs.m == m, "MATRIX SIZE MISMATCH");
            Matrix<T> res(n, m, 0);
            for (unsigned i = 0; i < n; i++) {
                for (unsigned j = 0; j < m; j++) {
                    res(i, j) = this->data[i][j] - rhs;
                }
            }
            return res;
        }

        //matrix addition
        void operator+=(const Matrix<T>& rhs) {
            ASSERT(rhs.n == n || rhs.m == m, "MATRIX SIZE MISMATCH");
            for (unsigned i = 0; i < n; i++) {
                for (unsigned j = 0; j < m; j++) {
                    this->data[i][j] += rhs(i, j);
                }
            }
        }

        //scalar addition
        void operator+=(const T& rhs) {
            ASSERT(rhs.n == n || rhs.m == m, "MATRIX SIZE MISMATCH");
            for (unsigned i = 0; i < n; i++) {
                for (unsigned j = 0; j < m; j++) {
                    this->data[i][j] += rhs;
                }
            }
        }

        //matrix subtraction
        void operator-=(const Matrix<T>& rhs) {
            ASSERT(rhs.n == n || rhs.m == m, "MATRIX SIZE MISMATCH");
            for (unsigned i = 0; i < n; i++) {
                for (unsigned j = 0; j < m; j++) {
                    this->data[i][j] -= rhs(i, j);
                }
            }
        }

        //scalar subtraction
        void operator-=(const T& rhs) {
            ASSERT(rhs.n == n || rhs.m == m, "MATRIX SIZE MISMATCH");
            for (unsigned i = 0; i < n; i++) {
                for (unsigned j = 0; j < m; j++) {
                    this->data[i][j] -= rhs;
                }
            }
        }

        //scalar multiplication
        Matrix<T> operator*(const T& rhs) {
            Matrix<T> res(this->n, this->m, 0);

            for (unsigned i = 0; i < this->m; i++) {
                for (unsigned j = 0; j < this->n; j++) {
                    res(i, j) += this->data[i][j] * rhs;
                }
            }
            return res;
        }

        //scalar multiplication
        void operator*=(const T& rhs) {
            Matrix<T> res = (*this) * rhs;
            (*this) = res;
        }

        //matrix multiplication
        Matrix<T> operator*(const Matrix<T>& rhs) {
            ASSERT(this->m == rhs.n, "MATRIX SIZE MISMATCH");
            Matrix<T> res(this->n, rhs.m, 0);
            for (unsigned i = 0; i < rhs.m; i++) {
                for (unsigned j = 0; j < this->n; j++) {
                    for (unsigned k = 0; k < rhs.n; k++) {
                        res(i, j) += this->data[i][k] * rhs(k, j);
                    }
                }
            }
            return res;
        }

        //matrix multiplication
        void operator*=(const Matrix<T>& rhs) {
            Matrix<T> res = (*this) * rhs;
            (*this) = res;
        }




        /*
        * Matrix Inversion Using Cholesky Decomposition Based on https://arxiv.org/ftp/arxiv/papers/1111/1111.4144.pdf
        */

        static Matrix<T> invOfCD(const Matrix<T>& A) {
            auto R = cd(A);
            //back-substitution

            for (int i = (A.n - 1); i >= 0; i--) {
                R(i, i) = (1 / R(i, i));
                for (int j = i - 1; j >= 0; j--) {
                    T sum = 0;
                    for (int k = i; k > j; k--) {
                        sum -= R(j, k) * R(k, i);
                    }
                    R(j, i) = (sum / R(j, j));
                }
            }

            return R;
        }

        //Cholesky Decomposition to get upper triangular matrix R
        static Matrix<T> cd(const Matrix<T>& A) {
            ASSERT(A.n == A.m, "MATRIX NOT SQUARE");
            Matrix<T> R(A.n, A.n, 0);
            for (unsigned i = 0; i < R.n; i++) {
                for (unsigned j = i; j < R.n; j++) {
                    T sum = 0;

                    //get diagonal elemets
                    if (i == j) {
                        for (unsigned k = 0; k < i; k++) {
                            sum += (R(k, i) * R(k, i));
                        }
                        R(i, i) = sqrt(A(i, i) - sum);
                    }

                    //get upper triangular values
                    else {
                        for (unsigned k = 0; k < i; k++) {
                            sum += (R(k, i) * R(k, j));
                        }
                        R(i, j) = (A(i, j) - sum) / R(i, i);
                    }
                }
            }
            return R;
        }

        static Matrix<T> cdInv(const Matrix<T>& A) {
            Matrix<T> R(A.n);
            auto Acdi = invOfCD(A);
            auto Acdit = Acdi.transpose();
            R = Acdi * Acdit;
            return R;
        }





        /*
        PSEUDO TODO
        addition dun
        subtraction dun
        multiplication dun
        determinant
        multiplication dun
        transpose dun
        inverse dun
        adjoint
        cofactor


        */



    };

}
