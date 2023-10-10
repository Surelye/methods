; Алгоритм разложения чисел в цепную дробь

(defun div (dividend divider)
  (multiple-value-bind (quotient)
      (floor dividend divider) quotient))


(defun chained-fraction-printer (chained-frac ps qs)
  (let* ((len (length ps))
         (indent (+ 2 (max (length (write-to-string (car (last ps))))
                           (length (write-to-string (car (last qs)))) 3)))
         (all-lists (list (cons 'i (loop for i from -1 to (- len 2) collect i))
                          (cons 'q_i (cons () (cons () chained-frac)))
                          (cons 'P_i ps) (cons 'Q_i qs))))
    (format t "~va" 4 " ")
    (do ((i 0 (1+ i))) ((= (* (1+ len) (1+ indent)) i) (terpri)) (format t "_"))
    (dolist (lst all-lists)
      (format t "~va|" 4 " ")
      (do ((i 0 (1+ i))) ((= (1+ len) i))
        (format t (format nil "~va" indent (nth i lst))) (format t "|"))
      (format t "~%~va|" 4 " ")
      (do ((i 1 (1+ i))) ((= (* (1+ len) (1+ indent)) i) (format t "|~%"))
        (format t "_")))))


(defun get-chained-fraction (numer denom)
  (let ((chained-fraction) (p-prev 0) (q-prev 1)
        (p-cur 1) (q-cur 0) (q-i) (ps '(1 0)) (qs '(0 1)))
    (aux::while (not (zerop (mod numer denom)))
      (setq q-i (div numer denom))
      (psetq chained-fraction (cons q-i chained-fraction)
             numer denom denom (mod numer denom)
             p-cur (+ (* q-i p-cur) p-prev) p-prev p-cur
             q-cur (+ (* q-i q-cur) q-prev) q-prev q-cur)
      (setq ps (cons p-cur ps) qs (cons q-cur qs)))
    (setq q-i (div numer denom) ps (cons (+ (* q-i p-cur) p-prev) ps)
          qs (cons (+ (* q-i q-cur) q-prev) qs)
          chained-fraction (reverse (cons q-i chained-fraction)))
    (chained-fraction-printer chained-fraction (reverse ps) (reverse qs))
    chained-fraction))


(defun chained-fraction-handler ()
  (format t "~%Введите числитель и знаменатель дроби через пробел: ")
  (let ((input))
    (tagbody try-again
       (setq input (uiop:split-string (setq input (read-line)) :separator " ")
             input (mapcar #'parse-integer input))
       (when (zerop (cadr input))
         (format t "Знаменатель дроби не может быть равен нулю! Попробуйте снова: ")
         (go try-again)))
    (format t "~%Таблица числителей и знаменателей подходящих дробей: ~%")
    (get-chained-fraction (car input) (cadr input))))


; Алгоритмы приложения цепных дробей
; Решение линейных диофантовых уравнений ax + by = c

(defun diophant-printer (a b c x y)
  (format t "Решение диофантова уравнения ~ax - ~ay = ~a имеет вид:~%" a b c)
  (format t "    x = (-1)^k * cQ_{k-1} + bt = ~a + ~at;~%" x b)
  (format t "    y = (-1)^k * cP_{k-1} + at = ~a + ~at.~%" y a))


(defun diophant-machinerie (a b c)
  (format t "~%Таблица числителей и знаменателей подходящих дробей:~% ")
  (let* ((p-prev 0) (q-prev 1) (p-cur 1) (q-cur 0)
         (chained-fraction (get-chained-fraction a b))
         (len-chained (1- (length chained-fraction)))
         (q-i) (factor) (x) (y))
    (do ((i 0 (1+ i))) ((= len-chained i))
      (setq q-i (nth i chained-fraction))
      (psetq p-cur (+ (* q-i p-cur) p-prev) p-prev p-cur
             q-cur (+ (* q-i q-cur) q-prev) q-prev q-cur))
    (setq factor (expt -1 (logand (1+ len-chained) 1))
          x (* factor c q-cur)
          y (* factor c p-cur)) (terpri)
    (diophant-printer a b c x y)))


(defun diophant-handler ()
  (format t "~%Введите коэффициенты a, b, c диофантова уравнения через пробел: ")
  (let ((input))
    (tagbody try-again
       (setq input (uiop:split-string (setq input (read-line)) :separator " ")
             input (mapcar #'parse-integer input))
       (destructuring-bind (a b c) input
         (when (/= 1 (gcd a b))
           (format t "НОД(a, b) должен быть равен единице! Попробуйте ввести коэффициенты снова: ")
           (go try-again))
         (diophant-machinerie a b c)))))


; Решение линейных сравнений ax = b (mod m)

(defun solve-congruence-machinerie (a b m)
  (let ((gcd-a-m (gcd a m)) (chained-fraction)
        (len-chained) (x))
    (unless (zerop (mod b gcd-a-m))
      (return-from solve-congruence-machinerie))
    (setq a (/ a gcd-a-m) b (/ b gcd-a-m) m (/ m gcd-a-m)
          chained-fraction (get-chained-fraction m a)
          len-chained (1- (length chained-fraction)))
    (do ((i 0 (1+ i)) (p-cur 1 (+ (* (nth i chained-fraction) p-cur) p-prev))
         (p-prev 0 p-cur))
        ((= len-chained i)
         (progn (setq x (mod (* (expt -1 (logand len-chained 1)) p-cur b) m))
                (mapcar #'(lambda (mult) (+ x (* mult m)))
                        (loop for mult from 0 to (1- gcd-a-m) collect mult)))))))


(defun congruence-printer (res a b m)
  (if (= 1 b)
      (format t "~%Обратный к ~a по модулю ~a элемент x равен ~a.~%"
              a m (car res))
      (format t "~%Сравнение ~a * x = ~a (mod ~a) имеет решение:
    x = (~{~a~^, ~}).~%" a b m res)))


(defun congruence-handler ()
  (format t "~%Введите коэффициенты сравнения a, b и модуль m через пробел: ")
  (let ((input) (res))
    (setq input (uiop:split-string (setq input (read-line)) :separator " ")
          input (mapcar #'parse-integer input))
    (destructuring-bind (a b m) input
      (setq res (solve-congruence-machinerie a b m))
      (when (null res)
        (format t "~%Сравнение с заданными параметрами не может быть решено!~%")
        (return-from congruence-handler))
      (congruence-printer res a b m))))


; Алгоритм вычисления символа Лежандра

(defun is-square-residue (a p)
  (= 1 (aux::mod-expt a (ash (1- p) -1) p)))


(defun compute-legendre-machinerie (a p)
  (cond ((zerop (mod a p)) 0)
        ((is-square-residue a p) 1)
        (t -1)))


(defun compute-legendre (a p)
  (when (and (integerp a) (aux::miller-rabin p))
    (compute-legendre-machinerie a p)))


; Алгоритм вычисления символа Якоби

(defun compute-jacobi-machinerie (a b)
  (let ((r 1) (t-val) (c))
    (when (< a 0) (setq a (- a))
      (when (= 3 (mod b 4)) (setq r (- r))))
    (tagbody eliminate-evenness
       (setq t-val 0)
       (aux::while (zerop (logand a 1))
         (setq t-val (1+ t-val) a (ash a -1)))
       (when (= 1 (logand t-val 1))
         (when (or (= 3 (mod b 8)) (= 5 (mod b 8)))
           (setq r (- r))))
       (when (and (= 3 (mod a 4)) (= 3 (mod b 4)))
         (setq r (- r)))
       (setq c a a (mod b c) b c)
       (if (not (zerop a))
           (go eliminate-evenness)
           (return-from compute-jacobi-machinerie r)))))


(defun compute-jacobi (a b)
  (when (and (integerp a) (integerp b) (= 1 (logand b 1)) (> b 1))
    (if (= 1 (gcd a b))
        (compute-jacobi-machinerie a b) 0)))


(defun symbols-printer (a p legendre jacobi)
  (if (null legendre)
      (format t "~%    Поскольку p = ~a не является простым, символ Лежандра не может быть вычислен." p)
      (format t "~%    Символ Лежандра (~a/~a) равен ~a." a p legendre))
  (format t "~%    Символ Якоби    (~a/~a) равен ~a.~%" a p jacobi))


(defun symbols-handler ()
  (format t "~%Введите a и p для вычисления символов через пробел: ")
  (let ((input) (legendre) (jacobi))
    (tagbody try-again
       (setq input (uiop:split-string (setq input (read-line)) :separator " ")
             input (mapcar #'parse-integer input))
       (destructuring-bind (a p) input
         (when (evenp p)
           (format t "~%Параметр p не может быть чётным числом! Введите a и p снова: ")
           (go try-again))
         (setq jacobi (compute-jacobi a p)
               legendre (when (aux::miller-rabin p 20) (compute-legendre a p)))
         (symbols-printer a p legendre jacobi)))))


; Алгоритмы извлечения квадратного корня в Zp

(defun find-k-i (a-i q p)
  (do ((k 0 (1+ k))) ((= 1 (aux::mod-expt a-i (* (expt 2 k) q) p)) k)))


(defun get-inv (a p)
  (cadr (aux::ext-gcd a p)))


(defun seq-sqrt-Zp (a p)
  (when (/= 1 (compute-legendre a p)) (return-from seq-sqrt-Zp))
  (let ((b) (k-i -1) (k-is) (r-i) (m 0) (q (1- p))
        (a-prev a) (a-cur a) (pow))
    (aux::while (zerop (logand q 1)) (setq m (1+ m) q (ash q -1)))
    (aux::while (/= -1 (compute-legendre (setq b (random p)) p)))
    (aux::while (not (zerop k-i))
      (setq  k-i    (find-k-i a-cur q p) k-is (cons k-i k-is))
      (psetq a-cur  (mod (* a-cur (aux::mod-expt b (ash 1 (- m k-i)) p)) p)
             a-prev a-cur))
    (setq k-is (cdr k-is) r-i (aux::mod-expt a-prev (ash (1+ q) -1) p))
    (do ((i (length k-is) (1- i))) ((= 0 i) r-i)
      (setq pow (ash 1 (- m (car k-is) 1))
            r-i (mod (* r-i (get-inv (aux::mod-expt b pow p) p)) p)
            k-is (cdr k-is)))))


; Вероятностный алгоритм Чипполы

(defun complex-*-mod (f s sqr p)
  (let ((ffc (car f)) (fsc (cadr f))
        (sfc (car s)) (ssc (cadr s)))
    (mapcar #'(lambda (c) (mod c p))
            (list (+ (* ffc sfc) (* fsc ssc sqr))
                  (+ (* ffc ssc) (* fsc sfc))))))


(defun lst-mod-m (lst m)
  (mapcar #'(lambda (c) (mod c m)) lst))


(defun complex-^-mod (base power sqr m)
  (setq base (lst-mod-m base m))
  (do ((product '(1 0)))
      ((zerop power) (car product))
    (do () ((oddp power))
      (setq base (lst-mod-m (complex-*-mod base base sqr m) m)
            power (ash power -1)))
    (setq product (lst-mod-m (complex-*-mod product base sqr m) m)
          power (1- power))))


(defun cipolla (n p)
  (when (or (not (integerp n)) (not (integerp p)) (not (aux::miller-rabin p))
            (< p 3) (/= 1 (compute-legendre n p))) (return-from cipolla))
  (let ((a) (sqr))
    (aux::while (/= -1 (compute-legendre (- (* (setq a (random p)) a) n) p)))
    (setq sqr (- (* a a) n))
    (complex-^-mod (list a 1) (ash (1+ p) -1) sqr p)))


(defun roots-printer (a p seq cip)
  (when (or (null seq) (null cip))
    (format t "~%Поскольку ~a не является квадратичным вычетом по модулю ~a, корень не может быть найден.~%"
            a p) (return-from roots-printer))
   (format t "~%Корень, найденный по алгоритму Чипполы равен: ~a. Проверка: ~a * ~a (mod ~a) = ~a."
          cip cip cip p (mod (* cip cip) p))
   (format t "~%Корень, найденный с использованием арифметики поля равен: ~a. Проверка: ~a * ~a (mod ~a) = ~a.~%"
          seq seq seq p (mod (* seq seq) p)))


(defun roots-handler ()
  (format t "~%Введите число a, корень которого хотите найти, и модуль p через пробел: ")
  (let ((input) (seq-res) (cip-res))
    (setq input (uiop:split-string (setq input (read-line)) :separator " ")
          input (mapcar #'parse-integer input))
    (destructuring-bind (a p) input
      (when (not (aux::miller-rabin p 20))
        (format t "~%Число p должно быть простым!")
        (return-from roots-handler))
      (setq seq-res (seq-sqrt-zp a p)
            cip-res (cipolla a p))
      (roots-printer a p seq-res cip-res))))


; Функция вызова

(defun caller ()
  (let ((c -1))
    (aux::while t
      (format t "~%[1] -- Разложение чисел в цепную дробь;
[2] -- Решение линейных диофантовых уравнений ax + by = c;
[3] -- Нахождение обратного элемента в Z_m, решение линейных сравнений ax = b (mod m);
[4] -- Символы Лежандра и Якоби;
[5] -- Алгоритмы извлечения квадратного корня в кольце вычетов;
[0] -- Выход.~%")
      (format t "~%Ваш выбор: ")
      (aux::while (not (member (setq c (read)) '(1 2 3 4 5 0)))
        (format t "Некорректный ввод! Попробуйте снова: "))
      (cond ((= 1 c) (chained-fraction-handler))
            ((= 2 c) (diophant-handler))
            ((= 3 c) (congruence-handler))
            ((= 4 c) (symbols-handler))
            ((= 5 c) (roots-handler))
            (t (return-from caller))))))
