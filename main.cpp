#include <iostream>
#include "Chord.h"
#include "Scripts.h"
#include <string>
#include <vector>

using namespace std;
int main()
{
    // Выбираем параметры расчетов
    double tau = 0.05,
            h = 0.05,
            L = 4*M_PI,
            t0 = 0.0,
            T = 5.2,
            a = 1.0;
    int node = L/h + 1,
        test_index = 4;
    /**
     * @param test_index - выбор теста для решения
     * 1 - первый тест из лабораторной работы
     * 2- второй тест из лабороторной работы
     * 3 - тест моего варианта
     * 4 - для Мишани
     */
    struct ParametrsStruct parametrsStruct = {tau, h, L, t0, T, a, node, test_index};
    Scripts::clearAll();
    Chord chord(parametrsStruct);
    chord.run();
    Scripts::drawPics(test_index);
    return 0;
}