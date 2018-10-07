#include <iostream>
#include <math.h>
bool isSolution = false;
//checks if there exists a solution
double i1, i2, i3, f1, f2, s1, s2, val1, val2;
int xshift, yshift;
/*
 i1, i2, i3: times that the hydrophone recorded the signal
 f1, f2: x and y coordinated for 1st equation
 s1, s2: x and y coordinates for 2nd equation
 val1, val2: relative distance between the hydrophones
 val1 = s_sound*(i2-i1)
 val2 = s_sound*(i3-i1)
 notice that val1 and val2 are directed lengths and can be negative
 */
double coord[3][2];
//initial coordinates serves little purpose other than input
double s_sound;
//speed of sound causes equations to vary

double orgfunc(double a, double b, double x, double y){
    //this returns the numerical value of the function
    //mostly used to check for extraneous solutions
    return sqrt((x-a)*(x-a)+(y-b)*(y-b))-sqrt(x*x+y*y);
}

double func1(double a, double b, double k, double x){
    //original function solved for y with positive square root term
    return sqrt(k*k*(a*a+b*b-k*k)*(a*a+b*b-k*k-4*a*x+4*x*x)/(b*b-k*k)/(b*b-k*k))/2-(4*a*a*b+4*b*b*b-4*b*k*k-8*a*b*x)/(2*(4*k*k-4*b*b));
}

double deriv1(double a, double b, double k, double x){
    //d/dx(func1)
    return (k*k*(8*x-4*a)*(a*a+b*b-k*k))/(4*(b*b-k*k)*(b*b-k*k)*sqrt((k*k*(a*a+b*b-k*k)*(a*a-4*a*x+b*b-k*k+4*x*x))/(b*b-k*k)/(b*b-k*k)))+(4*a*b)/(4*k*k-4*b*b);
}


double func2(double a, double b, double k, double x){
    //original function solved for y with negative square root term
    return -sqrt(k*k*(a*a+b*b-k*k)*(a*a+b*b-k*k-4*a*x+4*x*x)/(b*b-k*k)/(b*b-k*k))/2-(4*a*a*b+4*b*b*b-4*b*k*k-8*a*b*x)/(2*(4*k*k-4*b*b));
}

double deriv2(double a, double b, double k, double x){
    //d/dx(func2)
    return (-k*k*(8*x-4*a)*(a*a+b*b-k*k))/(4*(b*b-k*k)*(b*b-k*k)*sqrt((k*k*(a*a+b*b-k*k)*(a*a-4*a*x+b*b-k*k+4*x*x))/(b*b-k*k)/(b*b-k*k)))+(4*a*b)/(4*k*k-4*b*b);
}

double absolval(double x){
    //returns the absolute value of x
    return x>0?x:-x;
}

void assign(){
    //assigns values after inputing
    //look above for variables
    if(i1<=i2 && i1<=i3){
        val1 = s_sound*(i2-i1);
        f1 = coord[1][0]-coord[0][0];
        f2 = coord[1][1]-coord[0][1];
        val2 = s_sound*(i3-i1);
        s1 = coord[2][0]-coord[0][0];
        s2 = coord[2][1]-coord[0][1];
        xshift = coord[0][0];
        yshift = coord[0][1];
        xshift = 0; yshift = 0;
    }
    else if (i2<=i1 && i2<=i3){
        val1 = s_sound*(i1-i2);
        f1 = coord[0][0]-coord[1][0];
        f2 = coord[0][1]-coord[1][1];
        val2 = s_sound*(i3-i2);
        s1 = coord[2][0]-coord[1][0];
        s2 = coord[2][1]-coord[1][1];
        xshift = coord[1][0];
        yshift = coord[1][1];
    }
    else if (i3<=i1 && i3<=i2){
        val1 = s_sound*(i1-i3);
        f1 = coord[0][0]-coord[2][0];
        f2 = coord[0][1]-coord[2][1];
        val2 = s_sound*(i2-i3);
        s1 = coord[1][0]-coord[2][0];
        s2 = coord[1][1]-coord[2][1];
        xshift = coord[2][0];
        yshift = coord[2][1];
    }
}

void input(){
    std::cout<<"What is the speed of sound in this medium?"<<std::endl;
    std::cin>>s_sound;
    std::cout<<"What are the coordinates of the first hydrophone?"<<std::endl;
    std::cin>>coord[0][0]>>coord[0][1];
    std::cout<<"What is the time?"<<std::endl;
    std::cin>>i1;
    std::cout<<"What are the coordinates of the second hydrophone?"<<std::endl;
    std::cin>>coord[1][0]>>coord[1][1];
    std::cout<<"What is the time?"<<std::endl;
    std::cin>>i2;
    std::cout<<"What are the coordinates of the third hydrophone?"<<std::endl;
    std::cin>>coord[2][0]>>coord[2][1];
    std::cout<<"What is the time?"<<std::endl;
    std::cin>>i3;
}

int main(int argc, const char * argv[]) {
    double xval[4] = {0, 0, 0, 0};
    //xval is the array for x_n values for the 4 different types of functions
    double repeat[4] = {0, 0, 0, 0};
    //stores previous values (if there are enough repeated values, there is a root)
    bool check[4] = {0, 0, 0, 0};
    //checks if xval[i] = NaN which means there is no solutions
    int count = 4;
    //number of functions that are not equal to NaN
    int repeatcount = 1;
    //number of times the same values have been repeated
    bool control = true;
    //a master control of the while loop, which is used to break it
    int increment = 0;
    //for more accuracy, we want a significantly large amount of iterations
    input(); assign();
    while(control || increment<(1<<12)){
        increment++;
        for(int k = 0; k<4; k++){
            if((xval[k]!=xval[k]) && check[k] == false) {
                //is xval[k] = NaN, then it is not equal to itself.
                count--;
                check[k] = true;
                if(count == 1){
                    control = false;
                    break;
                    //if only one variable left, control would be false
                    //although increment may still make the loop run
                }
            }
        }
        std::cout<<xval[0]<<" "<<xval[1]<<" "<<xval[2]<<" "<<xval[3]<<std::endl;
        xval[0] -= (func1(f1, f2, val1, xval[0])-func1(s1, s2, val2, xval[0]))/(deriv1(f1, f2, val1, xval[0])-deriv1(s1, s2, val2, xval[0]));
        xval[1] -= (func1(f1, f2, val1, xval[1])-func2(s1, s2, val2, xval[1]))/(deriv1(f1, f2, val1, xval[1])-deriv2(s1, s2, val2, xval[1]));
        xval[2] -= (func2(f1, f2, val1, xval[2])-func1(s1, s2, val2, xval[2]))/(deriv2(f1, f2, val1, xval[2])-deriv1(s1, s2, val2, xval[2]));
        xval[3] -= (func2(f1, f2, val1, xval[3])-func2(s1, s2, val2, xval[3]))/(deriv2(f1, f2, val1, xval[3])-deriv2(s1, s2, val2, xval[3]));
        //updates the xvalues according to the rule x_n+1 = x_n-f(x)/f'(x)
        bool checkloop = true;
        //boolean is true if all function values are repeated from last iteration to a certain accuracy
        for(int x = 0; x<4; x++){
            if(absolval(xval[x] - repeat[x]) > 0.0000001){
                //accuracy can be exchanged for speed, or vice versa
                checkloop = false;
                break;
            }
        }
        if(checkloop){
            repeatcount++;
            if(repeatcount>1000){
                //1000 repeats will break the loop
                //again, this value can be changed depending on accuracy needed
                break;
            }
        }
        else{
            repeatcount = 1;
            //updates the repeat values
            for(int x = 0; x<4; x++){
                repeat[x] = xval[x];
            }
        }
    }
    for(int x = 0; x<4; x++){
        //this means xval[x] != NaN
        if(!check[x]){
            /*each conditional checks for three things for the respective functions
             1. the y-value of the first equation is equal to the y-value of the second equation
             2. after plugging in the x and y values to the first the original equation, the right distance is calculated
             3. after plugging in the x and y values to the second the original equation, the right distance is calculated
             Note that 2 and 3 are solely there to check for outliers
             */
            if(x == 0 && absolval(func1(f1, f2, val1, xval[x]) - func1(s1, s2, val2, xval[x]))<0.000000001 &&
               absolval(val1 - orgfunc(f1, f2, xval[x], func1(f1, f2, val1, xval[x])))<0.00000001 &&
               absolval(val2 - orgfunc(s1, s2, xval[x], func1(f1, f2, val1, xval[x])))<0.00000001){
                isSolution = true;
                std::cout<<"("<<xval[x]+xshift<<", "<<func1(f1, f2, val1, xval[x])+yshift<<")"<<std::endl;
            }
            else if(x == 1 &&  absolval(func1(f1, f2, val1, xval[x]) - func2(s1, s2, val2, xval[x]))<0.000000001 &&
                    absolval(val1 - orgfunc(f1, f2, xval[x], func1(f1, f2, val1, xval[x])))<0.00000001 &&
                    absolval(val2 - orgfunc(s1, s2, xval[x], func1(f1, f2, val1, xval[x])))<0.00000001){
                isSolution = true;
                std::cout<<"("<<xval[x]+xshift<<", "<<func1(f1, f2, val1, xval[x])+yshift<<")"<<std::endl;
            }
            else if(x == 2 && absolval(func2(f1, f2, val1, xval[x]) - func1(s1, s2, val2, xval[x]))<0.0000000001 &&
                    absolval(val1 - orgfunc(f1, f2, xval[x], func2(f1, f2, val1, xval[x])))<0.00000001 &&
                    absolval(val2 - orgfunc(s1, s2, xval[x], func2(f1, f2, val1, xval[x])))<0.00000001){
                isSolution = true;
                std::cout<<"("<<xval[x]+xshift<<", "<<func2(f1, f2, val1, xval[x])+yshift<<")"<<std::endl;
            }
            else if(x == 3 && absolval(func2(f1, f2, val1, xval[x]) - func2(s1, s2, val2, xval[x]))<0.0000000001 &&
                    absolval(val1 - orgfunc(f1, f2, xval[x], func2(f1, f2, val1, xval[x])))<0.00000001 &&
                    absolval(val2 - orgfunc(s1, s2, xval[x], func2(f1, f2, val1, xval[x])))<0.00000001){
                isSolution = true;
                std::cout<<"("<<xval[x]+xshift<<", "<<func2(f1, f2, val1, xval[x])+yshift<<")"<<std::endl;
            }
        }
    }
    if(!isSolution){
        //if no root is found, which is only possible in unrealistic situations, "No Solution" is printed
        std::cout<<"No Solution"<<std::endl;
    }
}


// 1 -.5 0 0.000263474 .5 0 0.000707023 0 1 0
