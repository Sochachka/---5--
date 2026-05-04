import cmath
import math
import random
from fractions import Fraction
from typing import List, Union, Tuple

EPS = 1e-12
MAX_ITER = 1000

# ------------------------------------------------------------
# Ввод и базовые операции
# ------------------------------------------------------------
def parse_complex(s: str) -> complex:
    """Преобразует строку в комплексное число."""
    s = s.strip().lower().replace(' ', '')
    s = s.replace('i', 'j')
    if s == 'j':
        return 1j
    if s == '-j':
        return -1j
    try:
        return complex(s)
    except ValueError:
        raise ValueError(f"Неверный формат: {s}")

def input_complex(prompt: str) -> complex:
    while True:
        try:
            return parse_complex(input(prompt))
        except Exception as e:
            print(f"Ошибка: {e}\nПримеры: 5, -3.14, 3+4i, -2-5j, 7i")

# ------------------------------------------------------------
# Полином и производная
# ------------------------------------------------------------
def poly_val(coeffs: List[complex], x: complex) -> complex:
    """Вычисление полинома (коэффициенты от старшей степени к свободному)."""
    res = 0j
    for c in coeffs:
        res = res * x + c
    return res

def poly_der(coeffs: List[complex], x: complex) -> complex:
    """Производная."""
    n = len(coeffs)
    if n <= 1:
        return 0j
    res = 0j
    for i in range(n - 1):
        res = res * x + coeffs[i] * (n - i - 1)
    return res

def poly_scale_to_int(coeffs: List[complex]):
    """Приводит комплексные коэффициенты к целым числам, если возможно,
    возвращает масштабированные целые коэффициенты и множитель.
    Если коэффициенты не целые, возвращает исходные."""
    # Пытаемся представить каждое число как рациональное
    try:
        fracs = [Fraction(c.real) for c in coeffs]
        # проверим, что все мнимые части равны 0
        if any(abs(c.imag) > EPS for c in coeffs):
            return coeffs, 1
        # приведём к общему знаменателю
        lcm_den = 1
        for f in fracs:
            lcm_den = lcm_den * f.denominator // math.gcd(lcm_den, f.denominator)
        int_coeffs = [int(f * lcm_den) for f in fracs]
        return int_coeffs, lcm_den
    except Exception:
        return coeffs, 1

# ------------------------------------------------------------
# Рациональные корни (для вещественных целых коэффициентов)
# ------------------------------------------------------------
def rational_roots(coeffs: List[complex]) -> List[float]:
    """Возвращает список рациональных корней для полинома с целыми коэффициентами."""
    # Проверяем, что все коэффициенты вещественные (и целые)
    int_coeffs, scale = poly_scale_to_int(coeffs)
    if isinstance(int_coeffs[0], int) and all(isinstance(c, int) for c in int_coeffs):
        # Нормализуем, чтобы старший коэффициент был положительным
        a0 = int_coeffs[-1]  # свободный член
        an = int_coeffs[0]   # старший
        if an == 0:
            return []
        # Делители свободного члена
        def divisors(n):
            n = abs(n)
            divs = set()
            for i in range(1, int(math.sqrt(n)) + 1):
                if n % i == 0:
                    divs.add(i)
                    divs.add(n // i)
            return sorted(divs)
        p_divs = divisors(a0)
        q_divs = divisors(an)
        candidates = []
        for p in p_divs:
            for q in q_divs:
                candidates.append(p / q)
                candidates.append(-p / q)
        # Убираем дубликаты
        candidates = sorted(set(candidates))
        roots = []
        for r in candidates:
            if abs(poly_val(coeffs, complex(r, 0))) < 1e-8:
                roots.append(r)
        return roots
    return []

# ------------------------------------------------------------
# Численный метод Ньютона
# ------------------------------------------------------------
def newton(coeffs: List[complex], x0: complex) -> complex:
    """Метод Ньютона для комплексного корня."""
    x = x0
    for _ in range(MAX_ITER):
        fx = poly_val(coeffs, x)
        if abs(fx) < EPS:
            return x
        dfx = poly_der(coeffs, x)
        if abs(dfx) < EPS:
            break
        x_new = x - fx / dfx
        if abs(x_new - x) < EPS:
            return x_new
        x = x_new
    return x0

def find_roots_numerical(coeffs: List[complex]) -> List[complex]:
    """Находит все корни численно (понижая степень)."""
    roots = []
    # копируем коэффициенты
    cur = coeffs[:]
    # удаляем ведущие нули
    while len(cur) > 1 and abs(cur[0]) < EPS:
        cur.pop(0)
    if len(cur) == 1:
        return []
    # порядок
    n = len(cur) - 1
    # начальные приближения: на сетке + случайные
    guesses = []
    # вещественные на сетке
    for x in range(-3, 4):
        guesses.append(complex(x, 0))
    # комплексные
    for re in range(-2, 3):
        for im in range(-2, 3):
            if re != 0 or im != 0:
                guesses.append(complex(re, im))
    # случайные
    for _ in range(20):
        guesses.append(complex(random.uniform(-2, 2), random.uniform(-2, 2)))

    while len(roots) < n and cur:
        best = None
        # перебираем приближения
        for g in guesses:
            r = newton(cur, g)
            if is_root_valid(cur, r):
                if best is None or abs(poly_val(cur, r)) < abs(poly_val(cur, best)):
                    best = r
        if best is None:
            # последняя попытка: взять случайное
            best = newton(cur, complex(random.uniform(-2, 2), random.uniform(-2, 2)))
        roots.append(best)
        # делим полином
        cur = poly_div(cur, best)
        while len(cur) > 1 and abs(cur[0]) < EPS:
            cur.pop(0)
    return roots

def poly_div(coeffs: List[complex], root: complex) -> List[complex]:
    """Деление полинома на (x - root) (схема Горнера)."""
    n = len(coeffs)
    if n <= 1:
        return []
    new = [0j] * (n - 1)
    new[0] = coeffs[0]
    for i in range(1, n - 1):
        new[i] = new[i-1] * root + coeffs[i]
    return new

def is_root_valid(coeffs: List[complex], r: complex, eps=1e-8) -> bool:
    return abs(poly_val(coeffs, r)) < eps

# ------------------------------------------------------------
# Точные решения для квадратного уравнения (с радикалами)
# ------------------------------------------------------------
def solve_quadratic_exact(a, b, c):
    """Возвращает корни квадратного уравнения a x^2 + b x + c = 0.
    Если a, b, c целые, пытается представить в радикалах.
    Возвращает список корней в виде строк (если точно) или чисел."""
    # преобразуем в рациональные числа
    try:
        a_frac = Fraction(a.real) if abs(a.imag) < EPS else None
        b_frac = Fraction(b.real) if abs(b.imag) < EPS else None
        c_frac = Fraction(c.real) if abs(c.imag) < EPS else None
    except:
        a_frac = b_frac = c_frac = None
    if a_frac and b_frac and c_frac:
        # пытаемся вычислить дискриминант
        D = b_frac*b_frac - 4*a_frac*c_frac
        # если D отрицательный, корни комплексные, выводим через sqrt
        # попробуем извлечь квадрат из D как рациональное число
        sqrt_D = None
        # если D точный квадрат
        if D >= 0:
            # проверяем, является ли D полным квадратом рационального числа
            # числитель и знаменатель
            num = D.numerator
            den = D.denominator
            # ищем квадратные корни
            sqrt_num = math.isqrt(num) if num >= 0 else None
            sqrt_den = math.isqrt(den) if den >= 0 else None
            if sqrt_num and sqrt_den and sqrt_num*sqrt_num == num and sqrt_den*sqrt_den == den:
                sqrt_D = Fraction(sqrt_num, sqrt_den)
        if sqrt_D is not None:
            # корни рациональны
            r1 = (-b_frac + sqrt_D) / (2*a_frac)
            r2 = (-b_frac - sqrt_D) / (2*a_frac)
            # преобразуем в строки с дробями
            return [format_root_str(r1), format_root_str(r2)]
        else:
            # выводим с радикалами
            return [f"({-b_frac} ± √({D}))/({2*a_frac})"]
    # иначе численно
    disc = b*b - 4*a*c
    sqrt_disc = cmath.sqrt(disc)
    r1 = (-b + sqrt_disc) / (2*a)
    r2 = (-b - sqrt_disc) / (2*a)
    return [r1, r2]

def format_root_str(r):
    """Форматирует рациональное число в строку."""
    if r.denominator == 1:
        return str(r.numerator)
    return f"{r.numerator}/{r.denominator}"

# ------------------------------------------------------------
# Основная функция решения
# ------------------------------------------------------------
def solve_equation(coeffs: List[complex]):
    """Выводит корни и проверку."""
    # Определяем степень
    degree = len(coeffs) - 1
    # Если все коэффициенты равны нулю
    all_zero = all(abs(c) < EPS for c in coeffs)
    if all_zero:
        print("\nУравнение: 0 = 0")
        print("Ответ: любое число")
        return

    while degree >= 0 and abs(coeffs[0]) < EPS:
        coeffs.pop(0)
        degree -= 1
    if degree == 0:
        if abs(coeffs[0]) < EPS:
            print("\nУравнение: 0 = 0")
            print("Ответ: любое число")
        else:
            print("\nУравнение: ненулевая константа")
            print("Ответ: нет решений")
        return

    # Вывод уравнения
    print("\nУравнение:", end=' ')
    terms = []
    for i, c in enumerate(coeffs):
        if abs(c) > EPS:
            if i == 0:
                term = f"{format_num(c)} x^{degree}"
            elif i == degree-1:
                term = f"{format_num(c)} x"
            elif i == degree:
                term = f"{format_num(c)}"
            else:
                term = f"{format_num(c)} x^{degree-i}"
            terms.append(term)
    print(" + ".join(terms).replace("+ -", "- "), "= 0")

    # Ищем рациональные корни
    rationals = rational_roots(coeffs)
    roots = []
    # Если есть рациональные, используем их как точные
    if rationals:
        for r in rationals:
            roots.append(complex(r, 0))
        # понижаем степень
        cur = coeffs[:]
        for r in rationals:
            cur = poly_div(cur, complex(r, 0))
            while len(cur) > 1 and abs(cur[0]) < EPS:
                cur.pop(0)
        # остальные корни численно
        if len(cur) > 1:
            numeric = find_roots_numerical(cur)
            roots.extend(numeric)
    else:
        # нет рациональных – решаем численно
        roots = find_roots_numerical(coeffs)

    # Убираем дубликаты (близкие корни)
    unique = []
    for r in roots:
        if not any(abs(r - u) < 1e-6 for u in unique):
            unique.append(r)

    # Вывод корней (без проверки)
    print("\n=== Корни ===")
    for i, r in enumerate(unique, 1):
        if abs(r.imag) < 1e-10:
            # вещественный
            val = r.real
            # пробуем представить как дробь
            frac = to_fraction(val)
            if frac:
                print(f"x{i} = {frac}")
            else:
                # убираем лишние нули
                print(f"x{i} = {val:.5f}".rstrip('0').rstrip('.'))
        else:
            # комплексный
            print(f"x{i} = {format_complex(r)}")

def format_num(z: complex) -> str:
    """Форматирует комплексное число для отображения в уравнении."""
    if abs(z.imag) < EPS:
        return format_float(z.real)
    else:
        return format_complex(z)

def format_float(x: float) -> str:
    """Форматирует вещественное число (если возможно как дробь)."""
    frac = to_fraction(x)
    if frac:
        return frac
    return f"{x:.5f}".rstrip('0').rstrip('.')

def to_fraction(x: float, max_den=100) -> Union[str, None]:
    """Пытается представить число как простую дробь."""
    try:
        frac = Fraction(x).limit_denominator(max_den)
        if abs(float(frac) - x) < 1e-8:
            if frac.denominator == 1:
                return str(frac.numerator)
            return f"{frac.numerator}/{frac.denominator}"
    except:
        pass
    return None

def format_complex(z: complex) -> str:
    real = z.real
    imag = z.imag
    if abs(real) < EPS:
        real = 0.0
    if abs(imag) < EPS:
        imag = 0.0
    real_str = to_fraction(real) or f"{real:.5f}".rstrip('0').rstrip('.')
    imag_abs = abs(imag)
    imag_str = to_fraction(imag_abs) or f"{imag_abs:.5f}".rstrip('0').rstrip('.')
    sign = '+' if imag >= 0 else '-'
    if real == 0:
        return f"{sign}{imag_str}i" if sign == '+' else f"-{imag_str}i"
    else:
        return f"{real_str} {sign} {imag_str}i"

# ------------------------------------------------------------
# Главная функция
# ------------------------------------------------------------
def main():
    random.seed()
    print("Введите коэффициенты уравнения:")
    print("a x⁵ + b x⁴ + c x³ + d x² + e x + f = 0\n")
    a = input_complex("a (x⁵): ")
    b = input_complex("b (x⁴): ")
    c = input_complex("c (x³): ")
    d = input_complex("d (x²): ")
    e = input_complex("e (x): ")
    f = input_complex("f (const): ")

    # Формируем список коэффициентов от старшей степени к свободному
    coeffs = [a, b, c, d, e, f]

    solve_equation(coeffs)

if __name__ == "__main__":
    main()