//
//  MyVector.h
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/25/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#ifndef hawkesInCpp_MyVector_h
#define hawkesInCpp_MyVector_h
#include <vector>
#include "Matrix.h"
#include "math.h"

using namespace std;

template<class T>
class MyVector: public Matrix<T>{
public:
    MyVector(int Row_Size) : Matrix<T>(Row_Size,1){
        length = Row_Size;
    }
    
    MyVector(Matrix<T> mat) :Matrix<T>(mat.Get_RowSize(),1){
        for(int i = 0; i < mat.Get_RowSize(); i++){
            Set_Element(i,mat.Get_Element(i,0));
        }
        length = mat.Get_RowSize();
    }
    
    void Set_Length(int Row_Size){
        Matrix<T>::Set_Size(Row_Size,1);
        length = Row_Size;
    }
    void Set_Element(int Row, T val){
        int col = 0;
        Matrix<T>::Set_Element(Row,col,val);
    }
    
    T Get_Element(int Row){
        int col = 0;
        T ret =  Matrix<T>::Get_Element(Row,col);
        return ret;
    }
    
    int GetLength(){
        return length;
    }
    
    T &operator() (int Row){
        int col = 0;
        T ret = Matrix<T>::Get_Element(Row,col);
        return ret;
    }
    
    MyVector<T> expo(){
        MyVector<T> ret(length);
        for(int i=0; i < length; i++){
            ret.Set_Element(i,exp(Get_Element(i)));
        }
        return ret;
    }
    
    MyVector<T> loga(){
        MyVector<T> ret(length);
        for(int i=0; i < length; i++){
            if(Get_Element(i) > 0)
                ret.Set_Element(i,log(Get_Element(i)));
            else{
                ret.Set_Element(i,0); cout << "something wrong with lambda.cond.int > 0" << endl;
            }
        }
        return ret;
    }
    
    T sum(){
        T ret = 0;
        for(int i=0; i < length; i++){
            ret = ret + (double) this->Get_Element(i);
        }
        return ret;
    }
    
    MyVector<bool> operator<=(const T val){
        MyVector<bool> temp(length);int count=0;
        for(int i=0; i < length; i++){
            if(Get_Element(i) <= val){
                temp.Set_Element(i,true); count++;
            }
            else{
                temp.Set_Element(i, false);
            }
        }
         return temp;
    }
    MyVector<bool> operator<(const T val){
        MyVector<bool> temp(length);int count=0;
        for(int i=0; i < length; i++){
            if(Get_Element(i) < val){
                temp.Set_Element(i,true); count++;
            }
            else{
                temp.Set_Element(i, false);
            }
        }
        return temp;
    }
    
    MyVector<T> operator[](MyVector<bool> indices){
        if(length != indices.GetLength()){
            cout << "length of vector and boolean indices are not same." << endl;
            exit(0);
        }
        
        int true_count = 0;
        for(int i=0; i< length; i++){
            true_count+= (indices(i));
        }
        MyVector<T> ret(true_count); int count = 0;
        for(int i=0; i < length; i++){
            if(indices(i)){
                ret.Set_Element(count,Get_Element(i));
                count++;
            }
        }
        return ret;
    }
    
    MyVector<T> AnotB(MyVector<T> B){
        vector<T> temp;
        for(int i = 0; i < GetLength(); i++){
            if(!B.AcontainsX(Get_Element(i)))
               temp.push_back(Get_Element(i));
        }
        MyVector<T> ret((int) temp.size());
        for(int j = 0; j < temp.size(); j++){
            ret.Set_Element(j,temp[j]);
        }
        return ret;
    }
    
    bool AcontainsX(T X){
        bool ret = false;
        for(int i=0; i < length; i++){
            if(Get_Element(i) == X)
                ret = true;
        }
        return ret;
    }
    
    int length;
    
};



#endif
