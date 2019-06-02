//
// Created by Коля on 2019-05-14.
//
#include <iostream>
#ifndef LAB3_SCRIPTS_H
#define LAB3_SCRIPTS_H


class Scripts {
public:
    static void clearAll(){
        system("../Results/ClearAll");
    }
    static void drawPics(int test_index){
        switch(test_index){
            case(1):
                system("../Results/Test1");
                break;
            case(2):
                system("../Results/Test2");
                break;
            case(3):
                system("../Results/Test3");
                break;
            case(4):
                system("../Results/Test4");
                break;
            case(5):
                system("../Results/Test5");
                break;
            case(6):
                system("../Results/Test6");
                break;
        }

    }
};


#endif //LAB3_SCRIPTS_H
