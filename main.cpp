#include <iostream>
#include "Chord.h"
#include "Scripts.h"
#include <string>
#include <vector>

using namespace std;
int main()
{
    // Выбираем параметры расчетов
    double h = 0.025,
            tau = 0.025,
            L = 1.0,
            t0 = 0.0,
            T = 1.2,
            a = 1.0;
    int node = L/h + 1,
        test_index = 1;
    /**
     * @param test_index - выбор теста для решения
     * 1 - первый тест из лабораторной работы
     * 2- второй тест из лабороторной работы
     * 3 - тест моего варианта
     */
    struct ParametrsStruct parametrsStruct = {tau, h, L, t0, T, a, node, test_index};
    Scripts::clearAll();
    Chord chord(parametrsStruct);
    chord.run();
    Scripts::drawPics(test_index);
    return 0;
}