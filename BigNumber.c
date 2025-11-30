
#include "BigNumber.h"
#include <stdio.h>
#include <time.h>

// Нормализация (убираем ведущие нули)
void bi_normalize(bigint *a) {
    // Пока длина > 1 и старшее слово равно 0, уменьшаем длину
    while (a->len > 1 && a->limb[a->len - 1] == 0)
        a->len--; 
}

// Обнуление большого числа
void bi_zero(bigint *a) {
    // Заполняем массив limb нулями
    memset(a->limb, 0, sizeof(a->limb));
    a->len = 1;
}

// Запись 32-битного числа в bigint 
void bi_from_uint(bigint *a, uint32_t v) {
    bi_zero(a); // обнуляем 
    a->limb[0] = v; // записываем значение в младшее слово
    a->len = 1; // длина = 1 слово
    bi_normalize(a); // убираем ведущие нули
}

// Копирование большого числа
void bi_copy(bigint *dst, const bigint *src) {
    // Копируем только занятые слова
    memcpy(dst->limb, src->limb, sizeof(uint32_t) * src->len);
    dst->len = src->len; // копируем длину
}

// Сравнение больших чисел
// Если a < b: -1, a > b: 1, a = b: 0
int bi_cmp(const bigint *a, const bigint *b) {
    // Сравниваем длины
    if (a->len < b->len) return -1;
    if (a->len > b->len) return 1;

    // Если длины равны, сравниваем слова, начиная со старшего
    for (int i = a->len - 1; i >= 0; --i) {
        if (a->limb[i] < b->limb[i]) return -1;
        if (a->limb[i] > b->limb[i]) return 1;
    }
    // Числа равны
    return 0;
}

// Проверка большого числа на равенство нулю
int bi_is_zero(const bigint *a) {
    // Одно слово, равное 0
    return (a->len == 1 && a->limb[0] == 0);
}

// c = a + b
void bi_add(bigint *c, const bigint *a, const bigint *b) {
    uint64_t carry = 0; // перенос
    int max; // максимальное количество слов, которые надо обработать

    // max = max(a->len, b->len)
    if (a->len > b->len) {
        max = a->len;
    } else {
        max = b->len;
    }

    // Идём по словам от младшего к старшему
    int i = 0;
    while (i < max || carry != 0) {

        uint64_t av = 0;    // текущее слово числа a
        uint64_t bv = 0;    // текущее слово числа b

        // Если слово существует в a — читаем его
        if (i < a->len) {
            av = a->limb[i];
        }

        // Если слово существует в b — читаем его
        if (i < b->len) {
            bv = b->limb[i];
        }

        // Складываем два слова и перенос
        uint64_t sum = av + bv + carry;

        // В результат записываем младшие 32 бита суммы
        c->limb[i] = (uint32_t)sum;

        // Новый перенос = старшие 32 бита суммы
        carry = sum >> 32; // сдвиг вправо

        i++;
    }

    // если перенос 0, длина = max
    // если перенос 1, длина = max + 1
    if (carry != 0) {
        c->len = max + 1;
        c->limb[c->len - 1] = (uint32_t)carry;   // добавляем перенос как старшее слово
    } else {
        c->len = max;
    }

    // Убираем ведущие нули
    bi_normalize(c);
}


// c = a + uint32
void bi_add_uint(bigint *c, const bigint *a, uint32_t v) {
    bigint tmp; // временная переменная
    bi_copy(&tmp, a);

    // Складываем младшее слово с v
    uint64_t sum = (uint64_t)tmp.limb[0] + v;
    tmp.limb[0] = (uint32_t)sum; // младшие 32 бита суммы
    uint64_t carry = sum >> 32; // перенос

    int i = 1; // начинаем со второго слова
    // Пока есть перенос и не вышли за длину
    while (carry && i < tmp.len) {
        uint64_t cur = (uint64_t)tmp.limb[i] + carry;
        tmp.limb[i] = (uint32_t)cur; 
        carry = cur >> 32;
        i++;
    }
    // Если перенос остался и есть место в массиве, добавляем новое слово
    if (carry && i < MAX_LIMBS) {
        tmp.limb[i] = (uint32_t)carry;
        tmp.len = i + 1;
    }

    bi_normalize(&tmp); // убираем ведущие нули
    bi_copy(c, &tmp);   // копируем результат в с
}

// c = a - b (a >= b)
void bi_sub(bigint *c, const bigint *a, const bigint *b) {
    int64_t carry = 0; // заём
    int max = a->len; // длина результата

    for (int i = 0; i < max; ++i) {
        int64_t av = a->limb[i]; // текущее слово а

        int64_t bv; // текущее слово b
        if (i < b->len) {
            bv = b->limb[i];   // берем i-е слово из b
        } else {
            bv = 0;            // если слова нет, считаем его равным нулю
        }

        // Вычисляем разность с учётом займа
        int64_t diff = av - bv - carry;
        if (diff < 0) { // если ушли в -
            diff += ((int64_t)1 << 32); // берём в долг у следующего слова
            carry = 1; // заняли
        } else {
            carry = 0;
        }
        c->limb[i] = (uint32_t)diff; // пишем младшие 32 бита результата
    }
    c->len = max;
    bi_normalize(c); // убираем ведущие нули
}

// Сдвиг влево на 1 бит
void bi_shl1(bigint *a) {
    uint32_t carry = 0; // перенос

    for (int i = 0; i < a->len; ++i) {
    // Собираем новое 64-битное значение: старое слово сдвигаем на 1 и добавляем перенос
        uint64_t v = ((uint64_t)a->limb[i] << 1) | carry;
        a->limb[i] = (uint32_t)v; // младшие 32 бита — новое слово
        carry = (uint32_t)(v >> 32); // старшие биты — новый перенос
    }
    // Если после прохода остался перенос — добавляем новое слово, если есть место
    if (carry && a->len < MAX_LIMBS) {
        a->limb[a->len++] = carry;
    }
}

// Сдвиг вправо на 1 бит
void bi_shr1(bigint *a) {
    uint32_t carry = 0; // перенос

    for (int i = a->len - 1; i >= 0; --i) {
    // Собираем 64-бит: перенос занимает старшие 32 бита, текущее слово — младшие
        uint64_t v = ((uint64_t)carry << 32) | a->limb[i];
        a->limb[i] = (uint32_t)(v >> 1); // делим на 2: берём верхние 32 бита после сдвига
        carry = (uint32_t)(v & 1); // младший бит - перенос
        if (i == 0) break; // чтобы не уйти в -1
    }
    bi_normalize(a); // убираем ведущие нули
}

// c = a * b 
void bi_mul(bigint *c, const bigint *a, const bigint *b) {
    uint64_t tmp[2 * MAX_LIMBS]; // результат
    memset(tmp, 0, sizeof(tmp)); // обнуляем результат

    // Проходим по всем словам числа a
    for (int i = 0; i < a->len; ++i) {
        uint64_t carry = 0;  // перенос

        // Умножаем a->limb[i] на каждое слово b
        int j = 0;
        while (j < b->len || carry != 0) {

            uint64_t bv = 0; // текущее слово из b 

            // Если j меньше длины b, берём b->limb[j]
            if (j < b->len) {
                bv = b->limb[j];
            } else {
                bv = 0;
            }

            // cur = то, что уже было в tmp[i + j]
            //     + произведение текущих слов a и b
            //     + перенос из предыдущего шага
            uint64_t cur = tmp[i + j] +
                           (uint64_t)a->limb[i] * (uint64_t)bv +
                           carry;

            // В tmp[i + j] записываем только младшие 32 бита результата
            tmp[i + j] = (uint32_t)cur;

            // В carry сохраняем старшие 32 бита (перенос для следующего разряда)
            carry = cur >> 32;

            j++; // Переходим к следующему слову b
        }
    }

    // длина результата
    int res_len = a->len + b->len;
    // Ограничиваем длину, чтобы не выйти за пределы массива
    if (res_len > MAX_LIMBS) res_len = MAX_LIMBS;

    // Перекладываем результат из tmp в c
    for (int i = 0; i < res_len; ++i) {
        c->limb[i] = (uint32_t)tmp[i];
    }
    c->len = res_len;
    bi_normalize(c); // убираем ведущие нули
}

// умножение на 32-битовое слово: c = a * w
void bi_mul_word(bigint *c, const bigint *a, uint32_t w) {
    uint64_t carry = 0; // перенос

    for (int i = 0; i < a->len; ++i) {
        // Умножаем каждое слово и добавляем перенос
        uint64_t cur = (uint64_t)a->limb[i] * w + carry;
        c->limb[i] = (uint32_t)cur; // записываем младшие 32 бита
        carry = cur >> 32; // старшие биты - перенос
    }
    c->len = a->len;
    if (carry && c->len < MAX_LIMBS) { // если есть перенос и место
        c->limb[c->len++] = (uint32_t)carry; // добавляем новое слово
    }
    bi_normalize(c); // убираем ведущие нули
}

// Количество бит в числе
int bi_bit_length(const bigint *a) {
    if (bi_is_zero(a)) return 0; // у нуля длина 0

    // Берём старшее слово
    uint32_t top = a->limb[a->len - 1];
    int bits = (a->len - 1) * 32; // битов в полных словах, кроме старшего
    int leading = 0; // битов занято в старшем слове

    // Считаем, сколько раз можем сдвинуть top вправо, пока он не станет 0
    while (top) {
        top >>= 1;
        leading++;
    }
    // Общая длина = полные слова + биты в старшем
    return bits + leading;
}

// Значение i-го бита
int bi_get_bit(const bigint *a, int i) {
    int limb = i / 32; // номер слова, где находится бит
    int off = i % 32; // смещение внутри слова
    if (limb >= a->len) return 0; // если вылезли за длину, бит = 0
    return (a->limb[limb] >> off) & 1U; // достаём нужный бит
}

// a mod m (m - uint32_t)
uint32_t bi_mod_word(const bigint *a, uint32_t m) {
    uint64_t rem = 0; // остаток

    // Идём от старшего слова к младшему
    for (int i = a->len - 1; i >= 0; --i) {
        uint64_t cur = (rem << 32) | a->limb[i]; // приписываем очередное слово
        rem = cur % m;  // новый остаток
    }
    return (uint32_t)rem;
}

// Деление БЧ на слово: q = a / d, *r - остаток
void bi_div_word(const bigint *a, uint32_t d, bigint *q, uint32_t *r) {
    uint64_t rem = 0; // остаток
    bigint res; // результат деления
    bi_zero(&res);
    res.len = a->len;

    // Идём от старшего слова к младшему
    for (int i = a->len - 1; i >= 0; --i) {
        uint64_t cur = (rem << 32) | a->limb[i]; // добавляем новое слово к остатку
        uint32_t qdigit = (uint32_t)(cur / d); // целая часть — это цифра результата
        rem = cur % d; // новый остаток
        res.limb[i] = qdigit; // сохраняем цифру в соответствующее слово
        // Чтобы i не ушёл в -1
        if (i == 0) break;
    }
    bi_normalize(&res); // убираем ведущие нули
    bi_copy(q, &res); // копируем результат в q

    // Если нужен остаток — возвращаем его через r
    if (r) *r = (uint32_t)rem;
}

// res = x mod m (остаток от деления большого числа x на большое число m)
void bi_mod(bigint *res, const bigint *x, const bigint *m) {
    bigint tmp; // остаток
    bi_zero(&tmp);
    int xbits = bi_bit_length(x); // количество бит у x

    // Идём по битам x от старшего к младшему
    for (int i = xbits - 1; i >= 0; --i) {
        // Остаток <<= 1 (умножаем на 2)
        bi_shl1(&tmp);

        if (bi_get_bit(x, i)) {  // если текущий бит x равен 1
            bigint one;
            bi_from_uint(&one, 1); // one = 1
            bi_add(&tmp, &tmp, &one); // tmp = tmp + 1
        } 
        // Если tmp >= m, вычитаем m, чтобы сохранить tmp в диапазоне [0, m-1]
        if (bi_cmp(&tmp, m) >= 0) {
            bigint r;
            bi_sub(&r, &tmp, m); // r = tmp - m
            bi_copy(&tmp, &r); // tmp = r
        }
    }
    bi_copy(res, &tmp); // копируем результат в res
}

// c = (a * b) mod m
void bi_mulmod(bigint *c, const bigint *a, const bigint *b, const bigint *m) {
    bigint t;
    bi_mul(&t, a, b); // t = a * b 
    bi_mod(c, &t, m); // c = t mod m
}

// res = base^exp mod mod
void bi_powmod(bigint *res, const bigint *base, const bigint *exp, const bigint *mod) {
    bigint result, a;
    bi_from_uint(&result, 1); // result = 1
    bi_mod(&a, base, mod); // a = base mod mod

    // Количество бит в экспоненте
    int bits = bi_bit_length(exp);
    // Быстрое возведение в степень
    for (int i = 0; i < bits; ++i) {
        // Если i-й бит экспоненты равен 1
        if (bi_get_bit(exp, i)) {
            bigint tmp;
            bi_mulmod(&tmp, &result, &a, mod); // result = (result * a) mod mod
            bi_copy(&result, &tmp);
        }
        bigint tmp2;
        bi_mulmod(&tmp2, &a, &a, mod); // a = (a * a) mod mod
        bi_copy(&a, &tmp2);
    }
    bi_copy(res, &result); // копируем результат в res
}