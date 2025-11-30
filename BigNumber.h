#ifndef BIGNUMBER_H
#define BIGNUMBER_H

#define M_PI 3.14159265358979323846
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <wchar.h>
#include <locale.h>
#include <stdint.h>
#include <time.h>


#define KEY_BITS   1024        // ДЛИНА КЛЮЧА. ДЛЯ 65 536 БИТ ПОСТАВЬТЕ 65536
#define LIMB_BITS  32

// bigint может хранить числа до примерно 2*KEY_BITS бит
#define MAX_LIMBS  ((2 * KEY_BITS + LIMB_BITS - 1) / LIMB_BITS + 2)

// ===== BIGINT ТИП =====
typedef struct {
    uint32_t limb[MAX_LIMBS];  // младшее слово в limb[0]
    int len;                   // число используемых слов
} bigint;

// ===== RSA =====

typedef struct {
    bigint n;
    bigint e;
    bigint d;
} rsa_key;

// ===== БАЗОВЫЕ ОПЕРАЦИИ BIGINT =====

void bi_normalize(bigint *a);

void bi_zero(bigint *a);

void bi_from_uint(bigint *a, uint32_t v);

void bi_copy(bigint *dst, const bigint *src);

int bi_cmp(const bigint *a, const bigint *b);

int bi_is_zero(const bigint *a);

// c = a + b
void bi_add(bigint *c, const bigint *a, const bigint *b);
// c = a + uint32
void bi_add_uint(bigint *c, const bigint *a, uint32_t v);

// c = a - b (предполагается a>=b)
void bi_sub(bigint *c, const bigint *a, const bigint *b);

// сдвиг влево на 1 бит
void bi_shl1(bigint *a);

// сдвиг вправо на 1 бит
void bi_shr1(bigint *a);

// c = a * b (полное произведение, до ~2*KEY_BITS бит)
void bi_mul(bigint *c, const bigint *a, const bigint *b);

// умножение на 32-битовое слово: c = a * w
void bi_mul_word(bigint *c, const bigint *a, uint32_t w);

// количество бит в числе
int bi_bit_length(const bigint *a);

// бит i (0..)
int bi_get_bit(const bigint *a, int i);

// a mod m (m - uint32_t)
uint32_t bi_mod_word(const bigint *a, uint32_t m);
// деление на слово: q = a / d, *r = остаток
void bi_div_word(const bigint *a, uint32_t d, bigint *q, uint32_t *r);

// res = x mod m (медленный битовый алгоритм)
void bi_mod(bigint *res, const bigint *x, const bigint *m);

// c = (a * b) mod m
void bi_mulmod(bigint *c, const bigint *a, const bigint *b, const bigint *m);

// res = base^exp mod mod
void bi_powmod(bigint *res, const bigint *base, const bigint *exp, const bigint *mod);



// ===== HEX ПРЕОБРАЗОВАНИЕ =====

int hex_val(char c);

void bi_from_hex(bigint *a, const char *hex);

void bi_print_hex(const bigint *a);
// ===== СЛУЧАЙНЫЕ ЧИСЛА (НЕ криптостойко, но без внешних библиотек) =====

int get_random_bytes(uint8_t *buf, size_t len);

// Случайный big-int заданной длины в битах
int bi_random_bits(bigint *a, int bits);

// ===== МОДУЛЬНАЯ ОБРАТНАЯ d = e^{-1} mod phi (e - 32-битовое) =====

// расширенный gcd для 32-битных
uint32_t inv_mod_uint32(uint32_t a, uint32_t m);

// d = e^{-1} mod phi, где e - uint32_t (65537)
// Ищем k: (1 + k*phi) делится на e. d = (1 + k*phi)/e
int modinv_small_e(bigint *d, uint32_t e, const bigint *phi);
// ===== ПРОВЕРКА ПРОСТОТЫ (Миллер–Рабин) =====

int is_even(const bigint *n);



int is_probable_prime(const bigint *n);
// Генерация простого числа заданной битовой длины
int generate_prime(bigint *p, int bits);


// Генерация ключа RSA (bits, e=65537)
int rsa_generate_key(rsa_key *key, int bits);

// Шифрование: c = m^e mod n
void rsa_encrypt(bigint *c, const bigint *m, const rsa_key *k);
// Расшифровка: m = c^d mod n
void rsa_decrypt(bigint *m, const bigint *c, const rsa_key *k);

void bi_from_string(bigint *a, const char *s);
#endif // BIGNUMBER_H