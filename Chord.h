//
// Created by Коля on 2019-05-13.
//
#include <string>
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#pragma once
#ifndef UNTITLED_CHORD_H
#define UNTITLED_CHORD_H

using namespace std;
struct ParametrsStruct{
    double tau, h, L, t0, T, a;
    int node, test_index;
};

class Chord {
public:
    Chord(ParametrsStruct parametrsStruct) {
        this->tau = parametrsStruct.tau;
        this->h = parametrsStruct.h;
        this->L = parametrsStruct.L;
        this->t0 = parametrsStruct.t0;
        this->T = parametrsStruct.T;
        this->a = parametrsStruct.a;
        this->node = parametrsStruct.node;
        this->test_index = parametrsStruct.test_index;
        for(int i = 0; i < node; i++){
            this->initial_speed.push_back(0.0);
            this-> initial_deviation.push_back(0.0);
            this->fxx.push_back(0.0);
            this->g.push_back(0.0);
        }
    }
    ~Chord(){};
/*
 ****************               Задаем начальные условия                ****************
 */
     void initialDeviation(){
         double temp, Fxx;
         switch(test_index) {
             case(1):{
                 for (int i = 0; i < node; i++) {
                     temp = sin(M_PI * i * h);
                     /*Fxx = -M_PI * M_PI * sin(M_PI * i * h);
                     Fxx = sin(M_PI * (i+1) * h) - 2*sin(M_PI * i * h) + sin(M_PI * (i-1) * h);
                     Fxx*=1/(h*h);*/
                     this->initial_deviation[i] = temp;
                     //this->fxx[i] = Fxx;
                 }
                 break;
             }
             case(2):{
                 for(int i = 0; i < node; i++){
                     temp = i*h*(1 - i*h);
                     //Fxx = (1 - 2*i*h);
                     //Fxx = (i+1)*h*(1 - (i+1)*h) - 2*i*h*(1 - i*h) + (i-1)*h*(1 - (i-1)*h);
                     this->initial_deviation[i] = temp;
                     //this->fxx[i] = Fxx;
                 }
                 break;
             }
             case(3):{
                 for(int i = 0; i < node; i++){
                     temp = 0.5*(1 + i*i*h*h);
                     this->initial_deviation[i] = temp;
                 }
                 break;
             }
             case(4):{
                 for(int i = 0; i < node; i++){
                     temp = 0.0;
                     this->initial_deviation[i] = temp;
                 }
                 break;
             }
             case(5):{
                 for(int i = 0; i < node; i++){
                     if (abs(-2 + i*h) <= 0.166666666) {
                         temp = 1.0;
                         this->initial_deviation[i] = temp;
                     }
                     else{
                         temp = 0.0;
                         this->initial_deviation[i] = temp;
                     }
                 }
                 break;
             }
             case(6):{
                 for(int i = 0; i < node; i++) {
                     temp = 0.0;
                     this->initial_deviation[i] = temp;
                 }
                 fileOutPutLayer(0, initial_deviation);
             }
         }
     }
    void FxxCalculations(){
         this->fxx[0] = 0;
         for(int i = 1; i < node-1; i++) {
             this->fxx[i] = (initial_deviation[i + 1] - 2 * initial_deviation[i] + initial_deviation[i - 1]) / (h * h);
         }
        this->fxx[node-1] = 0;
     }

     void initialSpeed(){
         double temp;
         switch (test_index) {
             default:{
                 for(int i = 0; i < node; i++){
                     temp = 0.0;
                     this->g[i] =  temp;
                 }break;
             }
             case(3):{
                 for(int i = 0; i < node; i++){
                     temp = i * h * sin(2 * i * h);
                     this->g[i] =  temp;
                 }break;
             }
             case(6):{
                 for(int i = 0; i < node; i++) {
                     if (abs(-1 + i * h) <= 0.25) {
                         temp = 1.0 - 5.0*abs(-1 + i * h);
                         this->initial_deviation[i] = temp;
                     } else {
                         temp = 0.0;
                         this->initial_deviation[i] = temp;
                            }
                 }
             }break;
         }
     }
     void initialDeviationSecond(){
         double temp = 0.0;
         for(int i = 0; i < node; i++){
             temp = initial_deviation[i] + tau * g[i] + tau*tau*a*a*fxx[i]/2.0;
             this->initial_speed[i] = temp;
         }
         fileOutPutLayer(1, initial_speed);
     }

        double psi(double t, double y_next){     //Функция для ГУ на правом конце
         switch (test_index){
             case(1):
                 return 0.0;
             case(2):
                 return 0.0;
             case(3):
                 return 1.0;
             case(4):
                 return -sin(2.0*t);
             case(5):
                 return y_next; //Для ГУ второго рода
             case(6):
                 return y_next;
         }

        }
        double fi(double t, double y_next) {      //Функция для ГУ на левом конце
            switch (test_index){
                case(1):
                    return 0.0;
                case(2):
                    return 0.0;
                case(3):
                    return 0.5 + 3.0*t;
                case(4):
                    return sin(2.0*t);
                case(5):
                    return y_next;
                case(6):
                    return y_next;
            }
        }

     void fileOutPutLayer(int layer_num, vector<double> y_next){
         auto layer_number = to_string(layer_num);
         ofstream LayerOutput("../Results/" + layer_number + ".txt");
         switch(test_index) {
             default: {
                 for (int i = 0; i < node; i++) {
                     LayerOutput << i * h << "\t\t" << y_next[i] << endl;
                 }
                 break;
             }
             case(5):{
                 for (int i = 0; i < node; i++) {
                     LayerOutput << -2.0 + i * h << "\t\t" << y_next[i] << endl;
                 }
                 break;
             }
             case(6):{
                 for (int i = 0; i < node; i++) {
                     LayerOutput << -1.0 + i * h << "\t\t" << y_next[i] << endl;
                 }
                 break;
             }
         }
     }
     void stabilityCheck(){
         if(a*tau/h <= 1.0){
             cout << "Число Курранта = " << a*tau/h <<", схема стабильна." << endl;
         } else cout << "Число Курранта = " << a*tau/h <<", схема нестабильна." << endl;

     }

/*
 * Записываем разностную схему
 */
     void chordCalcutions(){
         vector<double> y_previous(node), y_current(node), y_next(node);
         double t = t0 + 2.0*tau,
                courant_num = tau*tau*a*a/(h*h),
                temp = 0.0,
                temp1 = 0.0;
         int layer_num = 2;

          for(int i = 0; i < node; i++){
             y_previous[i] = initial_deviation[i];
          }
          for(int i = 1; i < node-1; i++){
            y_current[i] = initial_speed[i];
          }
            y_current[0] = fi(t-tau, y_current[0]);
            y_current[node-1] = psi(t-tau, y_current[node-2]);

         do {

             for(int i = 1; i < node - 1; i++){
                 temp = y_current[i+1] - 2 * y_current[i] + y_current[i-1];
                 temp1 = courant_num*temp + 2 * y_current[i] - y_previous[i];
                 y_next[i] = temp1;
             }
             y_next[0] = fi(t, y_next[1]);
             y_next[node - 1] = psi(t, y_next[node-2]);

             errorLabTest(y_next, t, layer_num);

             for(int i = 0; i < node; i++) {
                 y_previous[i] = y_current[i];
                 y_current[i] = y_next[i];
             }

             fileOutPutLayer(layer_num, y_next);

             layer_num++;
             t += tau;
         } while(t <= T);
     }

     void errorLabTest(vector<double> y_next, double t, double layer_num){
         double error, maxError = 0.0, x, coord;
         switch(test_index) {
             case (1): {
                 for (int i = 0; i < node; i++) {
                     error = abs(y_next[i] - sin(M_PI * i * h) * cos(M_PI * t));
                     if (error > maxError) {
                         maxError = error;
                         x = i*h;
                     }
                 }
                 if(abs(t - 0.3) < eps){
                     cout << "B момент времени t == " << t << " ошибка равна "<< maxError <<endl;
                     cout << "Максимум ошибки получен в точке " << x <<endl;
                 }
                 //cout << "Ошибка на " << layer_num << "-м слое равна " << maxError << endl;
                 break;
             }
             case (2): {
                 double  u, x;
                 int k = 0.0;
                 k = int(0.5 * sqrt(2 / (M_PI * M_PI * eps))-1.0);
                 //cout << "k = " << k << endl;
                 for (int j = 0; j < node; j++) {
                     x = j * h;
                     u = 0.0;
                     for (int i = 0; i <= k; i++) {
                         u += (1.0 / ((2 * i + 1) * (2 * i + 1) * (2 * i + 1))) * sin((2 * i + 1) * M_PI * x) *
                               cos((2 * i + 1) * M_PI * t);
                     }
                     u *= 8.0 / (M_PI * M_PI * M_PI);
                     error = abs(y_next[j] - u);
                     if (error > maxError) {
                         maxError = error;
                         coord = j*h;
                     }
                 }
                 if(abs(t - 0.3) < eps){
                     cout << "B момент времени t == " << t << " ошибка равна "<< maxError <<endl;
                     cout << "Максимум ошибки получен в точке " << coord <<endl;
                 }
                 //cout << "Ошибка на " << layer_num << "-м слое равна " << maxError << endl;
                 break;
             }
             case(3):
                 break;
             case(4):
                 break;
             case(5):
                 break;
             case(6):
                 break;
         }
     }
     void run(){
         stabilityCheck();
         initialDeviation();
         FxxCalculations();
         initialSpeed();
         initialDeviationSecond();
         chordCalcutions();
     }


private:
    double tau, h, L, t0, T, a;
     int node, test_index;
    vector<double>  initial_deviation;//первый слой, начальное отклонение
    vector<double>  initial_speed; //значения на втором слое
    vector<double>  fxx;//вторая производная начального распределения
    vector<double>  g; //начальная скорость струны
    const  double eps = pow(10,-10);
};


#endif //UNTITLED_CHORD_H