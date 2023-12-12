(defpackage #:discrete
  (:use #:cl)
  (:export #:gelfond-shanks
           #:rho-pollard
           #:index-method))


(in-package #:discrete)


(defun read-params ()
  (format t "~%Введите значения a, b, p через пробел: ")
  (mapcar #'parse-integer (uiop:split-string (read-line) :separator " ")))


(defun gelfond-shanks (a b p)
  (let (r a-pows a-1 a-1-pows-b temp k+ris)
    (setq r (1+ (floor (sqrt p)))
          a-pows (sort (mapcar #'(lambda (pow) (list pow (aux:mod-expt a pow p)))
                               (loop for pow from 1 to (1- r) collect pow))
                       #'(lambda (f s) (< (cadr f) (cadr s))))
          a-1 (mod (cadr (aux:ext-gcd (aux:mod-expt a r p) p)) p)
          a-1-pows-b (mapcar #'(lambda (pow) (mod (* (aux:mod-expt a-1 pow p) b) p))
                             (loop for pow from 0 below r collect pow)))
    (loop for j from 0 to (1- r)
          do (when (setq temp (car (member (nth j a-1-pows-b) a-pows :key #'cadr)))
               (setq k+ris (cons (+ (car temp) (* r j)) k+ris)))
          finally (return (if (null k+ris)
                              nil
                              (apply #'min k+ris))))))


(defun rho-pollard (a b p eps)
  (let ((m 1) sqrt-m k (1/3-p (floor (* 1/3 p))) (2/3-p (floor (* 2/3 p)))
        i s y-0 y-i y-2i alpha-i alpha-2i beta-i beta-2i d (temp (mod a p))
        x-coef)
    (loop
      (when (= 1 temp)
        (return))
      (setq m (1+ m)
            temp (mod (* temp a) p)))
    (setq sqrt-m (sqrt m)
          k (1+ (floor (sqrt (* 2 sqrt-m (log (/ eps)))))))
    (labels ((f (y)
               (mod (cond ((<= 1 y 1/3-p) (* b y))
                          ((<= 1/3-p y 2/3-p) (* y y))
                          (t (* a y))) p))
             (which-subset? (y)
               (cond ((<= 1 y 1/3-p) 1)
                     ((<= 1/3-p y 2/3-p) 2)
                     (t 3)))
             (next-alpha (y alpha)
               (let ((subs-# (which-subset? y)))
                 (mod (cond ((= 1 subs-#) alpha)
                            ((= 2 subs-#) (* 2 alpha))
                            (t (1+ alpha))) m)))
             (next-beta (y beta)
               (let ((subs-# (which-subset? y)))
                 (mod (cond ((= 1 subs-#) (1+ beta))
                            ((= 2 subs-#) (* 2 beta))
                            (t beta)) m)))
             (solve-congruence (a-cong b-cong m-cong d)
               (let (x)
                 (unless (zerop (mod b-cong d))
                   (return-from solve-congruence nil))
                 (setq a-cong (/ (mod a-cong m) d)
                       b-cong (/ (mod b-cong m) d)
                       m-cong (/ m-cong d)
                       x (mod (* b-cong (cadr (aux:ext-gcd a-cong m-cong))) m-cong))
                 (do ((j 0 (1+ j))) ((> j d) nil)
                   (when (= b (aux:mod-expt a x p))
                     (return-from solve-congruence x))
                   (setq x (+ x m-cong))))))
      (tagbody
       step-2
         (setq i 1
               s (random m)
               y-0 (aux:mod-expt a s p)
               y-i (f y-0)
               y-2i (f y-i)
               alpha-i (next-alpha y-0 s)
               alpha-2i (next-alpha y-i alpha-i)
               beta-i (next-beta y-0 0)
               beta-2i (next-beta y-i beta-i))
         (go step-4)
       step-3
         (setq i (1+ i)
               alpha-i (next-alpha y-i alpha-i)
               beta-i (next-beta y-i beta-i)
               y-i (f y-i)
               alpha-2i (next-alpha y-2i alpha-2i) beta-2i (next-beta y-2i beta-2i)
               y-2i (f y-2i)
               alpha-2i (next-alpha y-2i alpha-2i) beta-2i (next-beta y-2i beta-2i)
               y-2i (f y-2i))
         (go step-4)
       step-4
         (if (/= y-i y-2i)
             (if (< i k)
                 (go step-3)
                 (return-from rho-pollard nil))
             (go step-5))
       step-5
         (setq x-coef (- beta-i beta-2i)
               d (gcd x-coef m))
         (cond ((= 1 d) (return-from rho-pollard (mod (* (- alpha-2i alpha-i) (cadr (aux:ext-gcd x-coef m))) m)))
               ((<= 2 d sqrt-m) (if (null (setq temp (solve-congruence x-coef (- alpha-2i alpha-i) m d)))
                                    (return-from rho-pollard nil)
                                    (return-from rho-pollard temp)))
               (t (go step-2)))))))


(defun rho-pollard-handler (a b p)
  (let (epsilon)
    (format t "~%Введите значение эпсилон (ро-метод): ")
    (setq epsilon (read-from-string (read-line)))
    (rho-pollard a b p epsilon)))


(defun index-method (a b p base)
  (unless (aux:check-primitive a p)
    (format t "~%Элемент ~a не является образующим!~%" a)
    (return-from index-method nil))
  (let (x)
    (aux:write-to-file (list p a b base) "params")
    (uiop:run-program "python3 main.py")
    (handler-case (setq x (parse-integer (uiop:read-file-line "result")))
      (error ()
        (return-from index-method nil)))
    (uiop:run-program "rm params result") x))


(defun index-method-handler (a b p)
  (let (base)
    (format t "~%Введите значение B, по которому будет построена факторная база (индекс-метод): ")
    (setq base (parse-integer (read-line)))
    (index-method a b p base)))


(defun print-results (a b p res-gs res-rho res-im)
  (let (temp)
    (if (null res-gs)
        (format t "~%По методу Гельфонда-Шенкса значение дискретного логарифма найти не удалось.~%")
        (unless (= -1 res-gs)
          (format t "~%Согласно методу Гельфонда-Шенкса значение дискретного логарифма равно: ~a.
Проверка: ~a^~a = ~a (mod ~a), ~a = ~a: ~a.~%"
                  res-gs a res-gs (setq temp (aux:mod-expt a res-gs p)) p b temp (eql b temp))))
    (if (null res-rho)
        (format t "~%По rho-методу Полларда значение дискретного логарифма найти не удалось.~%")
        (unless (= -1 res-rho)
          (format t "~%Согласно rho-методу Полларда значение дискретного логарифма равно: ~a.
Проверка: ~a^~a = ~a (mod ~a), ~a = ~a: ~a.~%"
                  res-rho a res-rho (setq temp (aux:mod-expt a res-rho p)) p b temp (eql b temp))))
    (if (null res-im)
        (format t "~%По индекс-методу значение дискретного логарифма найти не удалось.~%")
        (unless (= -1 res-im)
          (format t "~%Согласно индекс-методу значение дискретного логарифма равно: ~a.
Проверка: ~a^~a = ~a (mod ~a), ~a = ~a: ~a.~%"
                  res-im a res-im (setq temp (aux:mod-expt a res-im p)) p b temp (eql b temp)))) t))

(defun caller ()
  (let ((choice -1) (res-gs -1) (res-rho -1) (res-im -1))
    (aux:while t
      (format t "~%[1] -- Метод Гельфонда-Шенкса вычисления дискретного логарифма в произвольной циклической группе;
[2] -- Rho-метод Полларда вычисления дискретного логарифма в конечной циклической группе;
[3] -- Индекс-метод дискретного логарифмирования в конечном простом поле;
[4] -- Все методы;
[0] -- Выход.~%")
      (format t "~%Ваш выбор: ")
      (aux:while (not (member (setq choice (read)) '(1 2 3 4 0)))
        (format t "Некорректный ввод! Попробуйте снова: "))
      (when (= 0 choice)
        (return-from caller t))
      (destructuring-bind (a b p) (read-params)
        (cond ((= 1 choice) (setq res-gs (gelfond-shanks a b p)))
              ((= 2 choice) (setq res-rho (rho-pollard-handler a b p)))
              ((= 3 choice) (setq res-im (index-method-handler a b p)))
              ((= 4 choice) (setq res-gs (gelfond-shanks a b p)
                                  res-rho (rho-pollard-handler a b p)
                                  res-im (index-method-handler a b p)))
              (t 'SKIP))
        (print-results a b p res-gs res-rho res-im)))))
