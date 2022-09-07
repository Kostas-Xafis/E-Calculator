#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/*
gcc calcAll.c -o calcAll
*/
struct Args {
    int threads;
    int expo;
    float base;
    int len_pow;
};

void *newline() {
    FILE* fp;
    fp = fopen("allPerf.txt", "a+");
    fseek(fp, 0, SEEK_END);
    fputs("\n", fp);
    fclose(fp);
}

int main(){    
    struct Args args[18];
    int inc = 0;
    //Million to 10 Million 
    for(int i = 1; i <= 3; i++){
        args[inc].threads = 4;
        args[inc].expo = 6;
        args[inc].base = pow(i*1061520, 1./6);
        args[inc].len_pow = 14;
        inc++;
    }
    for(int i = 4; i <= 10; i++){
        args[inc].threads = 8;
        args[inc].expo = 6;
        args[inc].base = pow(i*1061520, 1./6);
        args[inc].len_pow = 16;
        inc++;
    }

    //10 Million to 100 Million
    for(int i = 2; i <= 9; i++){
        args[inc].threads = 8;
        args[inc].expo = 7;
        args[inc].base = pow(i*10355293, 1./7);
        args[inc].len_pow = 17;
        inc++;
    }

    for (int i = 0; i < 100; i++) {
        printf("Digits: %d\n", (i+1) * ((int) pow(10, 5 + ((int) log10(10 + i)))));
        char * str = (char *) malloc(sizeof(char) * 30); 
        sprintf(str, "calc9 %d %d %f %d ", args[0].threads, args[0].expo, args[0].base, args[0].len_pow);
        // sprintf(str, "calc9 %d %d %f %d -s", args[i].threads, args[i].expo, args[i].base, args[i].len_pow);
        for(int j = 0; j < 10; j++) system(str);
        newline();
        printf("\n\n");
    }
    return 0;
}