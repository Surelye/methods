(defpackage #:aux
  (:use :cl)
  (:export #:while
           #:ext-gcd
           #:mod-expt
           #:miller-rabin))


(in-package #:aux)


(defmacro while (condition &body body)
  `(loop while ,condition
         do (progn ,@body)))


(defun ext-gcd (a b)
  (let ((s 0) (old-s 1) (r b) (old-r a)
        (quotient) (bezout-t))
    (while (not (zerop r))
           (setq quotient (floor old-r r))
           (psetq old-r r r (- old-r (* quotient r))
                  old-s s s (- old-s (* quotient s))))
    (if (zerop b) (setq bezout-t 0)
        (setq bezout-t (floor (- old-r (* old-s a)) b)))
    (list old-r old-s bezout-t)))


(defun mod-expt (base power divisor)
  (setq base (mod base divisor))
  (do ((product 1))
      ((zerop power) product)
    (do () ((oddp power))
      (setq base (mod (* base base) divisor)
            power (ash power -1)))
    (setq product (mod (* product base) divisor)
          power (1- power))))


(defun miller-rabin (num &optional (rounds 10))
  (when (or (= 2 num) (= 3 num)) (return-from miller-rabin t))
  (when (zerop (mod num 2)) (return-from miller-rabin))
  (let* ((n-pred (1- num)) (s 0) (t-val n-pred) (round-num 0) (a) (x))
    (macro::while (zerop (mod t-val 2)) (setq s (1+ s) t-val (ash t-val -1)))
    (tagbody next-iteration
       (macro::while (< round-num rounds)
         (setq a (+ 2 (random (- num 3)))
               x (mod-expt a t-val num))
         (when (or (= x 1) (= x n-pred))
           (setq round-num (1+ round-num))
           (go next-iteration))
         (dotimes (iter (1- s))
           (setq x (mod (* x x) num))
           (when (= 1 x) (return-from miller-rabin))
           (when (= n-pred x)
             (setq round-num (1+ round-num))
             (go next-iteration)))
         (return-from miller-rabin))
       (return-from miller-rabin t))))
