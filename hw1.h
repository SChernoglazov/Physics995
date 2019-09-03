#ifndef HOMEWORK1_HW1_H
#define HOMEWORK1_HW1_H

using namespace std;

double interpolation(const double x, const double* points, const int N, const double dx);
double interpolation(const double* points, const int N);
double smooth(const double x);
double jump(const double x);
int main();

#endif //HOMEWORK1_HW1_H
