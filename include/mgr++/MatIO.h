#pragma once
#include <stdexcept>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include <matio.h>


template <class Derived>
int read_mat(const char* matpath, const char* varname, Derived& matrix) {
    using namespace std;

    try {
        // Open mat file.
        mat_t* mat = Mat_Open(matpath, MAT_ACC_RDONLY);
        if (!mat) 
            throw logic_error(string("Failed to open file: ") + matpath);

        // Read variable.
        matvar_t* mat_var = 0;
        mat_var = Mat_VarRead(mat, varname);
        if (!mat_var) 
            throw logic_error(string("Variable not found: ") + varname);

        // Construct output matrix.
        unsigned x_size = mat_var->nbytes/mat_var->data_size;
        double* x = static_cast<double*>(mat_var->data);
        if (mat_var->rank != 2)
            throw logic_error(string("Unsupported dimensionality: ") + to_string(mat_var->rank));
        int rows = mat_var->dims[0];
        int cols = mat_var->dims[1];
        matrix.resize(rows, cols);
        Eigen::Map<Derived> map(x, rows, cols);
        matrix = map;

        // Close mat file.
        Mat_Close(mat);

        return 1;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 0;
    }
}

template <class Derived>
int write_mat(const char* matpath, const char* varname, Derived& matrix) {
    using namespace std;

    try {
        // Open mat file.
        mat_t* mat = Mat_Open(matpath, MAT_ACC_RDWR);
        if (!mat)
            throw logic_error(string("Failed to open file: ") + matpath);

        // Write variable.
        size_t dims[2] = {(size_t) matrix.rows(), (size_t) matrix.cols()};
        matvar_t* mat_var = Mat_VarCreate(varname, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, matrix.data(), 0);
        Mat_VarDelete(mat, varname);
        Mat_VarWrite(mat, mat_var, MAT_COMPRESSION_ZLIB);
        Mat_VarFree(mat_var);

        Mat_Close(mat);

        return 1;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 0;
    }    
}