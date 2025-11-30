
// gcc -o rsa main.c BigNumber.c rsa.c -lm
// ./rsa
#include "BigNumber.h"
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <time.h>
#include <string.h>


int main(void) {

    printf("Генерация RSA ключа %d бит...\n", KEY_BITS);

    srand((unsigned)time(NULL));

    rsa_key key;
    bi_zero(&key.n);
    bi_zero(&key.e);
    bi_zero(&key.d);

    if (!rsa_generate_key(&key, KEY_BITS)) {
        fprintf(stderr, "Ошибка генерации ключа\n");
        return 1;
    }

    printf("n (modulus) = 0x");
    bi_print_hex(&key.n);
    printf("\n\n");

    printf("e (public)  = 0x");
    bi_print_hex(&key.e);
    printf("\n\n");

    printf("d (private) = 0x");
    bi_print_hex(&key.d);
    printf("\n\n");

    // Тест: m < n
    bigint m, c, m2;
    bi_from_uint(&m, 123456789u);
    // bi_from_string(&m, "hello");

    printf("m (исходное) = 0x");
    bi_print_hex(&m);
    printf("\n");

    rsa_encrypt(&c, &m, &key);
    printf("c (шифр)     = 0x");
    bi_print_hex(&c);
    printf("\n");

    rsa_decrypt(&m2, &c, &key);
    printf("m (расш)     = 0x");
    bi_print_hex(&m2);
    printf("\n");

    return 0;
}
