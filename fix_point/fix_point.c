#include <stdio.h>

int main() {
    double num;
    printf("input 0.[DEC something] : Ex. 0.123456789\n");
    printf("Enter number : ");
    scanf("%lf", &num);

    // แปลงเป็น float <1
//    double num = input / 1000.0;

    printf("0.");
    for(int i = 0; i < 24; i++) {
        num *= 2;
        if(num >= 1.0) {
            printf("1");
            num -= 1.0;
        } else {
            printf("0");
        }
    }
    printf("\n");

    return 0;
}

