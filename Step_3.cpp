#include <fstream>
#include <iostream>
#import <math.h>
double a, b, c;
int sum;

int main() {
    std::ifstream harambe;
    harambe.open("test.in");
    std::ofstream lives;
    lives.open("test.out");
    for(int x = 0; x<10; x++){
        harambe>>a>>b>>c;
        sum+=(a*a+b*b+c*c>200*200);
        lives<<(int)sqrt(a*a+b*b+c*c)<<std::endl;
    }
    lives<<sum<<std::endl;
    harambe.close();
    lives.close();
}


