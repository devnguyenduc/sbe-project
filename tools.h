#ifndef __TOOLS_H__
#define __TOOLS_H__
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
using namespace std;

//
template <typename T>
void copy_array(T *&array_need_copy, T *array_copy, int size_array)
{
    for (int i = 0; i < size_array; i++)
    {
        array_need_copy[i] = array_copy[i];
    }
}

template <typename T>
void copy_array(T **&array_need_copy, T **array_copy, int rows, int cols)
{
    for (int j = 0; j < rows; j++)
    {
        for (int k = 0; k < cols; k++)
        {
            array_need_copy[j][k] = array_copy[j][k];
        }
    }
}
//
template <typename T>
void copy_and_delete(T *&array_need_copy, T *array_copy, int size_array){
    for (int i = 0; i < size_array; i++)
    {
        array_need_copy[i] = array_copy[i];
    }
    delete []array_copy;
}
//
template<typename T>
T sum_array(T *array_need_sum, int N){
    T sum = array_need_sum[0];
    for(int i = 1; i < N; i++){
        sum += array_need_sum[i];
    }
    return sum;
}
// default save array 1 dimension
template <typename T>
void save_csv(string filename, T data, string mode){
    cout << data << endl;
    ofstream myFile;
    myFile.open(filename);
    if(myFile.is_open()){
        myFile << to_string(abs(data)) << endl;
    }
    myFile.close();
    
}
template <typename T>
void save_csv(string filename, T *array_save, int rows, string mode)
{
    string split = ",";
    ofstream myFile;
    myFile.open(filename);

    if (myFile.is_open())
    {  
        myFile << "LABEL" << endl;
        for (int j = 0; j < rows; j++)
        {
            myFile << to_string(abs(array_save[j])) << endl;
        }
        
    }
    myFile.close();
}

// default save arrat 2 dimensions
template <typename T>
void save_csv(string filename, T **array_save, int rows, int cols, string mode)
{
    string split = ",";
    ofstream myFile;
    myFile.open(filename);
    if (myFile.is_open())
    {
        // Set the Label for data
        string labels = "LABEL";
        string write_label = "LABEL";
        for (int k = 0; k < cols; k++){
            write_label += split + labels + to_string(k + 1);
        }
        myFile << write_label << endl;
        
        for (int j = 0; j < rows; j++)
        {
            string write_data = "";
            for(int k = 0; k < cols; k++){
                write_data += split + to_string(abs(array_save[j][k]));
            }
            myFile << write_data << endl;
        }
    }
    myFile.close();
}

// template <typename T>
// void save_gnuplot(string filename, T **array_save, int rows, int cols, string mode){
//     string split = ",";
//     ofstream myFile;
//     myFile.open(filename);
//     if (myFile.is_open())
//     {
//         // // Set the Label for data
//         // string labels = "LABEL";
//         // string write_label = "LABEL";
//         // for (int k = 0; k < cols; k++){
//         //     write_label += split + labels + to_string(k + 1);
//         // }
//         // myFile << write_label << endl;
//         // Write data
//         for (int j = 0; j < rows; j++)
//         {
//             string write_data = "";
//             for(int k = 0; k < cols; k++){
//                 write_data += split + to_string(abs(array_save[j][k]));
//             }
//             myFile << write_data << endl;
//         }
//     }
//     myFile.close();
// }

#endif
