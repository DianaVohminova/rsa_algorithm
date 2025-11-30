#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "BigNumber.h"

//HEX ПРЕОБРАЗОВАНИЕ 

// Перевод hex в число 0-15
// разность ASCII кодов символов даёт цифровое знач-е символа
int hex_val(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
    if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
    return -1;
}

// Запись hex в bigint
void bi_from_hex(bigint *a, const char *hex) {
    bi_zero(a); // обнуляем bigint

    // Если строка начинается с "0x" или "0X", пропускаем эти два символа
    if (hex[0] == '0' && (hex[1] == 'x' || hex[1] == 'X'))
        hex += 2; // cдвигаем указатель hex на 2 символа вперёд

    size_t len = strlen(hex); // длина строки hex 
    // Проходим по каждому символу строки слева направо
    for (size_t i = 0; i < len; ++i) {
        int v = hex_val(hex[i]); // переводим символ hex в число
        if (v < 0) continue; // если символ не hex, пропускаем

        // Обрабатываем одну hex-цифру v
        // Алгоритм: a = a * 16 + v (сдвиг влево на 4 бита и прибавление)
        uint64_t carry = 0; // перенос при умножении
        // Идём по всем словам числа a
        for (int j = 0; j < a->len; ++j) {
            uint64_t cur = ((uint64_t)a->limb[j] << 4) + carry; // сдвигаем слово на 4 бита влево (умножение на 16) + перенос
            a->limb[j] = (uint32_t)cur;  // записываем младшие 32 бита обратно в слово
            carry = cur >> 32; // старшие биты - перенос
        }
        // Если после цикла остался перенос и есть место
        if (carry && a->len < MAX_LIMBS) {
            a->limb[a->len++] = (uint32_t)carry; // добавляем новое слово с переносом
        }

        // a += v
        uint64_t sum = (uint64_t)a->limb[0] + v; // складываем младшее слово с v 
        a->limb[0] = (uint32_t)sum;  // записываем младшие 32 бита суммы
        uint64_t c2 = sum >> 32; // старшие биты суммы - перенос
        int idx = 1; // Начинаем обрабатывать перенос со второго слова

        // Пока есть перенос c2 и не вышли за текущую длину числа
        while (c2 && idx < a->len) {
            uint64_t cur = (uint64_t)a->limb[idx] + c2; // прибавляем перенос к следующему слову
            a->limb[idx] = (uint32_t)cur; // пишем результат в слово
            c2 = cur >> 32; // обновляем перенос
            idx++; // переходим к следующему слову
        }
        // Если перенос остался и мы можем расширить число
        if (c2 && idx < MAX_LIMBS) {
            a->limb[idx] = (uint32_t)c2; // записываем перенос в новое слово
            a->len = idx + 1; // увеличиваем длину
        }
    }
    if (a->len == 0) a->len = 1; // если длина стала 0, выставляем её в 1
    bi_normalize(a); // убираем ведущие нули
}

// Печать bigint в формате hex
void bi_print_hex(const bigint *a) {
    // Если число = 0, печатаем 0
    if (bi_is_zero(a)) { 
        printf("0");
        return;
    }
    int i = a->len - 1; // Индекс старшего слова

     // Печатаем старшее слово без ведущих нулей (%X не добавляет ведущие нули)
    printf("%X", a->limb[i]);
    // Печатаем все остальные слова с ведущими нулями (по 8 hex-цифр на слово)
    for (i = i - 1; i >= 0; --i) {
        printf("%08X", a->limb[i]); // %08X - всегда 8 символов, дополняя нулями слева

    }
}

// СЛУЧАЙНЫЕ ЧИСЛА 

// Заполняем буфер buf длиной len псевдослучайными байтами
int get_random_bytes(uint8_t *buf, size_t len) {
    // Проходим по всем байтам буфера
    for (size_t i = 0; i < len; ++i) {
        // Берём младшие 8 бит от rand() и записываем в буфер
        buf[i] = (uint8_t)(rand() & 0xFF);
    }
    return 1;
}

// Случайный bigint заданной длины в битах
int bi_random_bits(bigint *a, int bits) {
    bi_zero(a); // обнуляем bigint
    int limbs = (bits + 31) / 32; // сколько 32-битных слов нужно для хранения bits бит
    // Если не помещаемся в наш bigint, возвращаем ошибку
    if (limbs > MAX_LIMBS) return 0;

    uint8_t tmp[4 * MAX_LIMBS]; // буфер байт под случайные данные
    size_t bytes = (size_t)limbs * 4; // нужно байт = количество слов * 4
    
    // Заполняем tmp случайными байтами; если не получилось, ошибка
    if (!get_random_bytes(tmp, bytes)) return 0;

    // Преобразуем каждые 4 байта в одно 32-битное слово limb
    for (int i = 0; i < limbs; ++i) {
        a->limb[i] = ((uint32_t)tmp[4*i]) | // младший байт
                     ((uint32_t)tmp[4*i+1] << 8) | // + следующий байт << 8
                     ((uint32_t)tmp[4*i+2] << 16) |// + следующий байт << 16
                     ((uint32_t)tmp[4*i+3] << 24); // + старший байт << 24
    }
    a->len = limbs;

    // гарантируем точную длину бит
    int extra_bits = limbs * 32 - bits;  // сколько "лишних" бит в старшем слове 
    if (extra_bits > 0) {
        uint32_t mask = 0xFFFFFFFFu >> extra_bits; // маска, которая оставляет только нужные биты
        a->limb[limbs - 1] &= mask; // обнуляем лишние старшие биты, чтобы длина была ровно bits
    }
    // ставим старший бит (чтобы точно bits)
    int highest_bit = (bits - 1) % 32; // номер старшего бита в старшем слове
    a->limb[limbs - 1] |= (1u << highest_bit); // принудительно ставим этот бит в 1
    
    // Устанавливаем младший бит (чтобы число было нечётное)
    a->limb[0] |= 1u;

    bi_normalize(a); // убираем ведущие нули
    return 1;
}

//  МОДУЛЬНАЯ ОБРАТНАЯ d = e^{-1} mod phi (e - 32-битовое) 

// Расширенный gcd для 32-битных
// Обратный элемент a по модулю m (т.е. a^{-1} mod m)
uint32_t inv_mod_uint32(uint32_t a, uint32_t m) {
    int64_t t = 0, newt = 1;
    int64_t r = m, newr = a;

    // Пока остаток не 0, продолжаем алгоритм
    while (newr != 0) {
        int64_t q = r / newr;  // q = целая часть от деления r / newr
        int64_t tmp = newt;
        newt = t - q * newt; // обновляем newt = t - q * newt
        t = tmp; // t = старое newt

        tmp = newr;
        newr = r - q * newr; // newr = r - q * newr
        r = tmp; // r = старое newr
    }
    if (r > 1) return 0; // если gcd(a, m) != 1 (r > 1), обратного элемента не существует
    if (t < 0) t += m; // приводим t к положительному представлению 
    return (uint32_t)t; // возвращаем t как a^{-1} mod m
}

// d = e^{-1} mod phi, где e - uint32_t (65537)
// Ищем k: (1 + k*phi) делится на e. d = (1 + k*phi)/e
int modinv_small_e(bigint *d, uint32_t e, const bigint *phi) {
    uint32_t p = bi_mod_word(phi, e); // p = phi mod e
    if (p == 0) return 0; // Если phi делится на e, обратной нет

    uint32_t inv_p = inv_mod_uint32(p, e);  // inv_p = p^{-1} mod e (обратный к p по модулю e)
    if (inv_p == 0) return 0;

    // k = -inv_p (mod e)
    uint32_t k = (e - inv_p) % e;

    // tmp = phi * k
    bigint tmp;
    bi_mul_word(&tmp, phi, k);

    // tmp2 = tmp + 1
    bigint tmp2;
    bi_add_uint(&tmp2, &tmp, 1);

    // d = tmp2 / e,  rem = tmp2 % e
    uint32_t rem;
    bi_div_word(&tmp2, e, d, &rem);
    if (rem != 0) return 0; // если остаток не 0, деление "не ровное", ошибка

    return 1;  // успех: d = e^{-1} mod phi
}

// ПРОВЕРКА ПРОСТОТЫ (Миллер–Рабин) 

// Набор оснований для теста Миллера–Рабина
uint32_t mr_bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};

// Проверка на чётность: возвращает 1, если число чётное, 0 если нечётное.
int is_even(const bigint *n) {
    // Смотрим младший бит младшего слова
    return !(n->limb[0] & 1u);
}

// Тест Миллера-Рабина
int is_probable_prime(const bigint *n) {
    bigint two, three;
    bi_from_uint(&two, 2);   // two = 2
    bi_from_uint(&three, 3); // three = 3

    if (bi_cmp(n, &two) < 0) return 0;      // 0 или 1 - не простые
    if (bi_cmp(n, &three) <= 0) return 1;   // 2 или 3 - простые

    if (is_even(n)) return 0; // если n чётное и > 2 - составное

    // Раскладываем n - 1 = d * 2^s, где d нечётное
    bigint n_minus_1, d;
    bigint one;
    bi_from_uint(&one, 1); // one = 1
    bi_sub(&n_minus_1, n, &one); // n_minus_1 = n - 1

    bi_copy(&d, &n_minus_1); // d = n - 1
    int s = 0; // степень двойки

    // Пока d чётное (младший бит = 0) и d != 0
    while (!bi_is_zero(&d) && ((d.limb[0] & 1u) == 0)) {
        bi_shr1(&d); // d = d / 2 (сдвиг вправо на 1 бит)
        s++; // увеличиваем степень двойки
    }

    // Перебираем все заданные основания mr_bases
    for (size_t i = 0; i < sizeof(mr_bases)/sizeof(mr_bases[0]); ++i) {
        uint32_t a32 = mr_bases[i];  // Берём  основание 
        bigint a;
        bi_from_uint(&a, a32); // a = a32 как bigint
        if (bi_cmp(&a, n) >= 0) continue;  // Если a >= n, тест с этим основанием бессмысленен

        bigint x;
        bi_powmod(&x, &a, &d, n); // x = a^d mod n

        // Если x == 1 или x == n-1, число проходит этот раунд теста
        if (bi_cmp(&x, &one) == 0 || bi_cmp(&x, &n_minus_1) == 0)
            continue;

        int cont = 0; 
        // Делаем s-1 шагов: x = x^2 mod n
        for (int r = 1; r < s; ++r) {
            bigint x2;
            bi_mulmod(&x2, &x, &x, n); // x2 = x * x mod n
            bi_copy(&x, &x2); // x = x2

             // Если x стало n-1, число возможно простое, выходим из цикла
            if (bi_cmp(&x, &n_minus_1) == 0) {
                cont = 1; // успешный раунд
                break;
            }
        }
        if (!cont) return 0; // если ни разу не получили n-1, число составное
    }
    return 1; // вероятно простое
}

// Генерация простого числа p заданной битовой длины bits
int generate_prime(bigint *p, int bits) {
    while (1) {
        // Генерируем случайное число нужной длины
        if (!bi_random_bits(p, bits)) return 0;
        p->limb[0] |= 1u;  // делаем p нечётным (ставим младший бит в 1)
        // Если p вероятно простое по тесту Миллера–Рабина
        if (is_probable_prime(p)) 
            return 1;
        // если нет — пробуем заново
    }
}


// RSA

// Генерация ключа RSA  длиной bits, e=65537
int rsa_generate_key(rsa_key *key, int bits) {
    int p_bits = bits / 2; // доля бит для простого числа p (половина)
    int q_bits = bits - p_bits;  // остаток бит для q (чтобы p*q ≈ 2^bits)

    uint32_t e_val = 65537u; // стандартное значение e
    bi_from_uint(&key->e, e_val);  // key->e = 65537
 
    bigint one;  // one = 1
    bi_from_uint(&one, 1);

    // Цикл, пока не получится корректная пара ключей
    while (1) {
        // Два простых числа p и q
        bigint p, q;
        bigint phi, p1, q1; 
        bi_zero(&p); bi_zero(&q);
        bi_zero(&phi); bi_zero(&p1); bi_zero(&q1);

        // 1. Генерируем два больших простых числа p и q
        if (!generate_prime(&p, p_bits)) continue;
        if (!generate_prime(&q, q_bits)) continue;

        // избегаем p == q, если они равны, генерируем заново
        if (bi_cmp(&p, &q) == 0) continue;

        // 2. n = p * q - модуль rsa
        bi_mul(&key->n, &p, &q); // key->n = p * q

        // 3. phi = (p-1)*(q-1)
        bi_sub(&p1, &p, &one); // p1 = p - 1
        bi_sub(&q1, &q, &one); // q1 = q - 1
        bi_mul(&phi, &p1, &q1); // phi = (p - 1) * (q - 1)

        // 4. проверяем нод(e, phi) == 1, чтобы e было взаимно простым с phi
        uint32_t phi_mod_e = bi_mod_word(&phi, e_val); // phi_mod_e = phi mod e
        uint32_t r0 = e_val, r1 = phi_mod_e; // r0 = e, r1 = phi mod e
        
        // Алгоритм Евклида для поиска нод(e, phi)
        while (r1 != 0) {
            uint32_t t = r0 % r1; // t = r0 mod r1
            r0 = r1; // r0 = r1
            r1 = t;  // r1 = t
        }
        // Если нод(e, phi) != 1
        if (r0 != 1) {
            // не взаимно просты — пробуем ещё раз
            continue;
        }

        // 5. d = e^{-1} mod phi - приватная экспонента
        if (!modinv_small_e(&key->d, e_val, &phi)) { // пытаемся вычислить d
            continue;
        }

        return 1; // ключ успешно сгенерирован
    }
}

// Шифрование: c = m^e mod n
// m — исходное сообщение, k — ключ с полем e и n
void rsa_encrypt(bigint *c, const bigint *m, const rsa_key *k) {
    bi_powmod(c, m, &k->e, &k->n);
}

// Расшифровка: m = c^d mod n
// c — шифртекст, m — результат, k — ключ с полем d и n
void rsa_decrypt(bigint *m, const bigint *c, const rsa_key *k) {
    bi_powmod(m, c, &k->d, &k->n);
}

void bi_from_string(bigint *a, const char *s) {
    bi_zero(a); // обнуляем bigint
    // Пока текущий символ строки не нулевой
    while (*s) {
        // Сдвигаем bigint влево на 8 бит (умножаем число на 256)
        uint64_t carry = 0; // перенос

        // Проходим по всем словам числа
        for (int i = 0; i < a->len; i++) {
            // Сдвигаем слово на 8 бит влево и добавляем перенос
            uint64_t cur = ((uint64_t)a->limb[i] << 8) | carry;
            // Сохраняем младшие 32 бита обратно в слово
            a->limb[i] = (uint32_t)cur;
            carry = cur >> 32; // старшие биты - перенос
        }
        // Если после сдвига остался перенос
        if (carry && a->len < MAX_LIMBS)
            a->limb[a->len++] = (uint32_t)carry; // добавляем новое слово

        // Добавляем в число текущий символ ===
        // a = a + ASCII-код символа
        bi_add_uint(a, a, (uint8_t)*s);
        // Переходим к следующему символу строки
        s++;
    }
}