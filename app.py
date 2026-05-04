import math
import cmath
import re
import json
import random
from fractions import Fraction
from flask import Flask, request, jsonify, render_template
from core import (
    EPS, parse_complex, poly_val, poly_der, poly_div,
    rational_roots, find_roots_numerical, solve_quadratic_exact,
    to_fraction, format_complex
)

app = Flask(__name__, template_folder='templates', static_folder='static')


def evaluate_expression(expr: str, x_value: complex) -> complex:
    """Вычисляет математическое выражение с поддержкой функций."""
    expr = expr.replace(' ', '')
    expr = expr.replace('π', str(math.pi))
    expr = expr.replace('pi', str(math.pi))
    
    def safe_sin(z):
        try:
            return cmath.sin(z)
        except:
            return math.sin(z.real) if isinstance(z, complex) else math.sin(z)
    
    def safe_cos(z):
        try:
            return cmath.cos(z)
        except:
            return math.cos(z.real) if isinstance(z, complex) else math.cos(z)
    
    def safe_tan(z):
        try:
            return cmath.tan(z)
        except:
            return math.tan(z.real) if isinstance(z, complex) else math.tan(z)
    
    def safe_cot(z):
        try:
            t = cmath.tan(z)
            if abs(t) < 1e-10:
                return 1e10
            return 1 / t
        except:
            t = math.tan(z.real) if isinstance(z, complex) else math.tan(z)
            if abs(t) < 1e-10:
                return 1e10
            return 1 / t
    
    def safe_ln(z):
        try:
            return cmath.log(z)
        except:
            return math.log(z.real) if isinstance(z, complex) and z.imag == 0 else cmath.log(z)
    
    def safe_log(z):
        try:
            return cmath.log10(z)
        except:
            return math.log10(z.real) if isinstance(z, complex) and z.imag == 0 else cmath.log10(z)
    
    def safe_exp(z):
        try:
            return cmath.exp(z)
        except:
            return math.exp(z.real) if isinstance(z, complex) else math.exp(z)
    
    def safe_sqrt(z):
        try:
            return cmath.sqrt(z)
        except:
            return math.sqrt(z.real) if isinstance(z, complex) and z.imag == 0 and z.real >= 0 else cmath.sqrt(z)
    
    def safe_abs(z):
        return abs(z)
    
    safe_dict = {
        'x': x_value,
        'sin': safe_sin,
        'cos': safe_cos,
        'tan': safe_tan,
        'tg': safe_tan,
        'cot': safe_cot,
        'ctg': safe_cot,
        'ln': safe_ln,
        'log': safe_log,
        'exp': safe_exp,
        'sqrt': safe_sqrt,
        'abs': safe_abs,
        'pi': math.pi,
        'e': math.e
    }
    
    try:
        result = eval(expr, {"__builtins__": {}}, safe_dict)
        return complex(result)
    except Exception as e:
        raise ValueError(f"Ошибка вычисления выражения: {e}")


def detect_periodic_solution(eq_str):
    """
    Определяет, является ли уравнение тригонометрическим,
    и возвращает формулу для корней, если возможно.
    """
    eq_str = eq_str.replace(' ', '').lower()
    
    # Разделяем на левую и правую часть
    if '=' in eq_str:
        left, right = eq_str.split('=', 1)
    else:
        return None
    
    # Проверяем простые тригонометрические уравнения
    # sin(x) = a
    sin_match = re.match(r'sin\(x\)\s*=\s*([-\d.]+|π/\d+|\d+π/?\d*)', eq_str)
    if sin_match:
        a_str = sin_match.group(1)
        try:
            if 'π' in a_str:
                a = eval(a_str.replace('π', str(math.pi)))
            else:
                a = float(a_str)
        except:
            a = None
        
        if a is not None and -1 <= a <= 1:
            # arcsin(a)
            angle = math.asin(a)
            # Проверяем, является ли угол красивым (π/6, π/4, π/3 и т.д.)
            angle_pi = angle / math.pi
            fractions = [(1/6, "π/6"), (1/4, "π/4"), (1/3, "π/3"), (1/2, "π/2"), 
                        (2/3, "2π/3"), (3/4, "3π/4"), (5/6, "5π/6")]
            
            angle_str = None
            for frac, name in fractions:
                if abs(angle_pi - frac) < 0.01:
                    angle_str = name
                    break
            
            if angle_str is None:
                angle_str = f"{angle:.4f}"
            
            return f"x = {angle_str} + 2π·k, k ∈ ℤ<br>x = π - {angle_str} + 2π·k, k ∈ ℤ"
    
    # cos(x) = a
    cos_match = re.match(r'cos\(x\)\s*=\s*([-\d.]+|π/\d+|\d+π/?\d*)', eq_str)
    if cos_match:
        a_str = cos_match.group(1)
        try:
            if 'π' in a_str:
                a = eval(a_str.replace('π', str(math.pi)))
            else:
                a = float(a_str)
        except:
            a = None
        
        if a is not None and -1 <= a <= 1:
            angle = math.acos(a)
            angle_pi = angle / math.pi
            fractions = [(1/6, "π/6"), (1/4, "π/4"), (1/3, "π/3"), (1/2, "π/2"),
                        (2/3, "2π/3"), (3/4, "3π/4"), (5/6, "5π/6")]
            
            angle_str = None
            for frac, name in fractions:
                if abs(angle_pi - frac) < 0.01:
                    angle_str = name
                    break
            
            if angle_str is None:
                angle_str = f"{angle:.4f}"
            
            return f"x = ±{angle_str} + 2π·k, k ∈ ℤ"
    
    # tan(x) = a
    tan_match = re.match(r'tan\(x\)\s*=\s*([-\d.]+|π/\d+|\d+π/?\d*)', eq_str)
    if tan_match:
        a_str = tan_match.group(1)
        try:
            if 'π' in a_str:
                a = eval(a_str.replace('π', str(math.pi)))
            else:
                a = float(a_str)
        except:
            a = None
        
        if a is not None:
            angle = math.atan(a)
            angle_pi = angle / math.pi
            fractions = [(1/4, "π/4"), (1/3, "π/3"), (1/6, "π/6"), (0, "0"), 
                        (1/2, "π/2"), (2/3, "2π/3"), (3/4, "3π/4")]
            
            angle_str = None
            for frac, name in fractions:
                if abs(angle_pi - frac) < 0.01:
                    angle_str = name
                    break
            
            if angle_str is None:
                angle_str = f"{angle:.4f}"
            
            return f"x = {angle_str} + π·k, k ∈ ℤ"
    
    return None


def parse_equation(eq_str: str):
    """Парсит уравнение. Если содержит функции, возвращает None и True."""
    eq_str = eq_str.replace(' ', '')
    
    function_patterns = ['sin', 'cos', 'tan', 'tg', 'cot', 'ctg', 'ln', 'log', 'exp', 'sqrt', 'abs']
    has_functions = any(pattern in eq_str.lower() for pattern in function_patterns)
    
    if has_functions:
        return None, True
    
    try:
        if '=' in eq_str:
            left, right = eq_str.split('=', 1)
        else:
            left, right = eq_str, '0'
        
        left_coeffs = parse_poly_terms(left)
        right_coeffs = parse_poly_terms(right)
        
        max_len = max(len(left_coeffs), len(right_coeffs))
        left_coeffs = [0j] * (max_len - len(left_coeffs)) + left_coeffs
        right_coeffs = [0j] * (max_len - len(right_coeffs)) + right_coeffs
        
        coeffs = [left_coeffs[i] - right_coeffs[i] for i in range(max_len)]
        
        while len(coeffs) > 1 and abs(coeffs[0]) < EPS:
            coeffs.pop(0)
        
        return coeffs, False
    except:
        return None, True


def parse_poly_terms(expr: str):
    """Разбирает строку с полиномом."""
    expr = expr.replace(' ', '')
    expr = expr.replace('π', str(math.pi))
    expr = expr.replace('i', 'j')

    if expr and expr[0] not in ('+', '-'):
        expr = '+' + expr

    terms = []
    start = 0
    for i, ch in enumerate(expr):
        if ch in ('+', '-') and i > 0:
            terms.append(expr[start:i])
            start = i
    terms.append(expr[start:])

    coeffs_dict = {}
    for term in terms:
        if not term:
            continue
        sign = 1
        if term[0] == '-':
            sign = -1
            term = term[1:]
        elif term[0] == '+':
            term = term[1:]

        if 'x' in term:
            if '^' in term:
                base, exp_str = term.split('^')
                exp = int(exp_str)
                if 'x' in base:
                    coeff_part = base.replace('x', '')
                else:
                    coeff_part = base
                if coeff_part == '':
                    coeff_part = '1'
            else:
                exp = 1
                if term == 'x':
                    coeff_part = '1'
                else:
                    coeff_part = term.replace('x', '')
                    if coeff_part == '':
                        coeff_part = '1'
        else:
            exp = 0
            coeff_part = term

        try:
            coeff_val = parse_complex(coeff_part) * sign
        except Exception as e:
            raise ValueError(f"Неверный коэффициент '{coeff_part}': {e}")

        coeffs_dict[exp] = coeffs_dict.get(exp, 0j) + coeff_val

    if not coeffs_dict:
        return [0j]
    max_deg = max(coeffs_dict.keys())
    coeffs = [coeffs_dict.get(d, 0j) for d in range(max_deg, -1, -1)]
    return coeffs


def find_roots_with_functions(equation_func, max_roots=6, search_range=(-8, 8)):
    """Находит корни уравнения f(x)=0."""
    from core import EPS, MAX_ITER
    
    roots = []
    
    guesses = []
    start, end = search_range
    step = 0.5
    x = start
    while x <= end:
        guesses.append(complex(x, 0))
        x += step
    
    for re in range(-3, 4):
        for im in range(-2, 3):
            guesses.append(complex(re, im))
    
    for _ in range(30):
        guesses.append(complex(random.uniform(start, end), random.uniform(-2, 2)))
    
    def numerical_derivative(f, x, h=1e-6):
        return (f(x + h) - f(x - h)) / (2 * h)
    
    def newton_with_func(f, x0):
        x = x0
        for _ in range(MAX_ITER):
            try:
                fx = f(x)
                if abs(fx) < EPS:
                    return x
                dfx = numerical_derivative(f, x)
                if abs(dfx) < EPS:
                    break
                x_new = x - fx / dfx
                if abs(x_new - x) < EPS:
                    return x_new
                x = x_new
            except:
                break
        return x0
    
    for g in guesses:
        try:
            r = newton_with_func(equation_func, g)
            if abs(equation_func(r)) < 1e-8:
                if abs(r.imag) < 1e-8 and (r.real < start - 1 or r.real > end + 1):
                    continue
                duplicate = False
                for u in roots:
                    if abs(r - u) < 1e-6:
                        duplicate = True
                        break
                if not duplicate:
                    roots.append(r)
                    if len(roots) >= max_roots:
                        break
        except:
            continue
    
    real_roots = [r for r in roots if abs(r.imag) < 1e-8]
    complex_roots = [r for r in roots if abs(r.imag) >= 1e-8]
    real_roots.sort(key=lambda r: r.real)
    
    return real_roots + complex_roots


def get_roots_with_functions(equation_func, eq_str):
    """Возвращает список строк с корнями для функционального уравнения."""
    # Сначала проверяем, есть ли красивая формула
    periodic_formula = detect_periodic_solution(eq_str)
    if periodic_formula:
        # Находим несколько примеров для наглядности
        roots = find_roots_with_functions(equation_func, max_roots=4, search_range=(-8, 8))
        examples = []
        for r in roots[:4]:
            if abs(r.imag) < 1e-10:
                examples.append(f"{r.real:.4f}")
        if examples:
            return [periodic_formula, f"Примеры: {', '.join(examples)}"]
        return [periodic_formula]
    
    # Если нет красивой формулы, ищем численно
    roots = find_roots_with_functions(equation_func, max_roots=6, search_range=(-8, 8))
    
    result = []
    for r in roots:
        if abs(r.imag) < 1e-10:
            val = r.real
            frac = to_fraction(val)
            if frac:
                result.append(frac)
            else:
                result.append(f"{val:.6f}".rstrip('0').rstrip('.'))
        else:
            result.append(format_complex(r))
    
    unique = []
    for r in result:
        if r not in unique:
            unique.append(r)
    
    if not unique:
        return ["Корни не найдены"]
    return unique


def get_roots(coeffs):
    """Возвращает список строк с корнями для полинома."""
    while len(coeffs) > 1 and abs(coeffs[0]) < EPS:
        coeffs.pop(0)

    degree = len(coeffs) - 1
    if degree <= 0:
        if degree == 0 and abs(coeffs[0]) < EPS:
            return ["Любое число"]
        else:
            return ["Нет решений"]

    rationals = rational_roots(coeffs)
    roots = []
    cur = coeffs[:]
    for r in rationals:
        roots.append(complex(r, 0))
        cur = poly_div(cur, complex(r, 0))
        while len(cur) > 1 and abs(cur[0]) < EPS:
            cur.pop(0)

    if len(cur) == 3:
        a, b, c = cur
        quad_roots = solve_quadratic_exact(a, b, c)
        for root in quad_roots:
            if isinstance(root, str):
                roots.append(root)
            else:
                if abs(root.imag) < 1e-10:
                    val = root.real
                    frac = to_fraction(val)
                    if frac:
                        roots.append(frac)
                    else:
                        roots.append(f"{val:.6f}".rstrip('0').rstrip('.'))
                else:
                    roots.append(format_complex(root))
    elif len(cur) == 2:
        a, b = cur
        root = -b / a
        if abs(root.imag) < 1e-10:
            val = root.real
            frac = to_fraction(val)
            if frac:
                roots.append(frac)
            else:
                roots.append(f"{val:.6f}".rstrip('0').rstrip('.'))
        else:
            roots.append(format_complex(root))
    elif len(cur) > 1:
        numeric = find_roots_numerical(cur)
        roots.extend(numeric)

    unique = []
    for r in roots:
        if isinstance(r, str):
            if r not in unique:
                unique.append(r)
        else:
            if not any(isinstance(u, str) is False and abs(r - u) < 1e-6 for u in unique):
                unique.append(r)

    result = []
    for r in unique:
        if isinstance(r, str):
            result.append(r)
        else:
            if abs(r.imag) < 1e-10:
                val = r.real
                frac = to_fraction(val)
                if frac:
                    result.append(frac)
                else:
                    result.append(f"{val:.6f}".rstrip('0').rstrip('.'))
            else:
                result.append(format_complex(r))
    
    if not result:
        return ["Корни не найдены"]
    return result


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/about')
def about():
    return render_template('about.html')


@app.route('/instructions')
def instructions():
    return render_template('instructions.html')


@app.route('/solve', methods=['POST'])
def solve():
    data = request.get_json()
    if not data or 'equation' not in data:
        return jsonify({'success': False, 'error': 'Не передано уравнение'})

    eq_str = data['equation'].strip()
    if not eq_str:
        return jsonify({'success': False, 'error': 'Пустое уравнение'})

    try:
        coeffs_or_none, use_numerical = parse_equation(eq_str)
        
        if use_numerical:
            if '=' in eq_str:
                left, right = eq_str.split('=', 1)
                def equation_func(x):
                    return evaluate_expression(left, x) - evaluate_expression(right, x)
            else:
                def equation_func(x):
                    return evaluate_expression(eq_str, x)
            
            roots = get_roots_with_functions(equation_func, eq_str)
            degree = "?"
        else:
            coeffs = coeffs_or_none
            degree = len(coeffs) - 1
            roots = get_roots(coeffs)
            
    except Exception as e:
        return jsonify({'success': False, 'error': f'Ошибка: {str(e)}'})

    if not roots:
        roots = ["Корни не найдены"]

    if len(roots) == 1 and roots[0] in ("Любое число", "Нет решений"):
        message = roots[0]
        roots = []
    else:
        message = None

    return jsonify({
        'success': True,
        'degree': degree,
        'roots': roots,
        'message': message
    })


if __name__ == '__main__':
    app.run(debug=True)