(defpackage #:primality
  (:use :cl))

(in-package #:primality)


; Тест Ферма проверки чисел на простоту

(defun fermat (n &optional (k 10))
  (let ((n-pred (1- n)) (a))
    (do ((i 0 (1+ i))) ((> i k) t)
      (aux::while (/= 1 (gcd (setq a (1+ (random n-pred))) n)))
      (when (/= 1 (aux::mod-expt a n-pred n))
        (return-from fermat)))))


; Тест Соловея-Штрассена

(defun solovay-strassen (n &optional (k 10))
  (cond ((= n 2) t)
        ((or (< n 2) (zerop (logand n 1))) nil)
        (t (let ((n-prev (1- n)) (jacobi) (a) (expt-res)
                 (power (ash (1- n) -1)))
             (do ((i 0 (1+ i))) ((> i k) t)
               (if (= 1 (gcd (setq a (1+ (random n-prev))) n))
                   (progn (setq jacobi (aux::compute-jacobi a n)
                                expt-res (aux::mod-expt a power n))
                          (unless (or (= 1 jacobi expt-res)
                                      (and (= -1 jacobi) (= n-prev expt-res)))
                            (return-from solovay-strassen)))
                   (return-from solovay-strassen)))))))


; Тест Миллера-Рабина

(defun miller-rabin (n &optional (k 10))
  (when (or (= 2 n) (= 3 n)) (return-from miller-rabin t))
  (when (or (< n 2) (= 0 (logand n 1))) (return-from miller-rabin))
  (let* ((n-pred (1- n)) (bound (- n-pred 2)) (t-val n-pred) (s 0) (round 0) (x))
    (aux::while (= 0 (logand t-val 1)) (setq s (1+ s) t-val (ash t-val -1)))
    (do () (nil)
      (tagbody next-iteration
         (when (= k round) (return-from miller-rabin t))
         (setq x (aux::mod-expt (+ 2 (random bound)) t-val n))
         (when (or (= 1 x) (= n-pred x))
           (incf round) (go next-iteration))
         (do ((iter 0 (1+ iter))) ((= iter (1- s)) (return-from miller-rabin))
           (setq x (mod (* x x) n))
           (when (= 1 x) (return-from miller-rabin))
           (when (= n-pred x)
             (incf round) (go next-iteration)))))))


; Функция вызова

(defun primality-printer (n fr sr tr)
  (format t "~%Результаты проверки числа ~a на простоту тестами:
    Ферма:             ~a;
    Соловея-Штрассена: ~a;
    Миллера-Рабина:    ~a.~%" n fr sr tr)
  (when (and fr (not sr) (not tr))
    (format t "~%Возможно, число ~a является числом Кармайкла.~%" n)))


(defun primality-handler ()
  (format t "~%Введите число для проверки и количество раундов проверки через пробел: ")
  (let ((input) (fermat) (sol-strass) (mill-rab))
    (tagbody try-again
       (setq input (uiop:split-string (setq input (read-line)) :separator " ")
             input (mapcar #'parse-integer input))
       (destructuring-bind (n k) input
         (when (or (< n 2) (< k 1))
           (format t "Некорректный ввод! Попробуйте снова: ")
           (go try-again))
         (setq fermat (fermat n k) sol-strass (solovay-strassen n k)
               mill-rab (miller-rabin n k))
         (primality-printer n fermat sol-strass mill-rab)))))


(defun caller ()
  (let ((c -1))
    (aux::while t
      (format t "~%[1] -- Проверка числа на простоту;
[0] -- Выход.~%")
      (format t "~%Ваш выбор: ")
      (aux::while (not (member (setq c (read)) '(1 0)))
        (format t "Некорректный ввод! Попробуйте снова: "))
      (cond ((= 1 c) (primality-handler))
            (t (return-from caller))))))
