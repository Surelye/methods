(defpackage #:aux
  (:use #:cl)
  (:export #:while
           #:write-to-file
           #:ext-gcd
           #:mod-expt
           #:miller-rabin
           #:rho-pollard
           #:check-primitive))


(in-package #:aux)


(defmacro while (condition &body body)
  `(loop while ,condition
         do (progn ,@body)))


(defun write-to-file (data filename)
  (with-open-file (out filename :direction :output :if-exists :supersede
                                :if-does-not-exist :create)
    (dolist (param data)
      (format out "~a " param))))


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
    (while (zerop (mod t-val 2)) (setq s (1+ s) t-val (ash t-val -1)))
    (tagbody next-iteration
       (while (< round-num rounds)
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


(defun rho-pollard-machinerie (n x-0 &optional (c 1))
  (when (miller-rabin n) (return-from rho-pollard-machinerie 'PRIME))
  (let ((mapping (lambda (x) (mod (+ c (* x x)) n)))
        (a x-0) (b x-0) (q))
    (tagbody map
       (setq a (funcall mapping a)
             b (funcall mapping (funcall mapping b))
             q (gcd (- a b) n))
       (cond ((and (< 1 q) (< q n)) (return-from rho-pollard-machinerie
                                      (list q (miller-rabin q))))
             ((= n q) (return-from rho-pollard-machinerie))
             (t (go map))))))


(defun rho-pollard (n x-0)
  (let ((c 1) (head) (factor) (factors))
    (while (zerop (logand n 1))
      (setq factors (cons 2 factors) n (ash n -1)))
    (setq x-0 (mod x-0 n))
    (while (/= 1 n)
      (setq factor (rho-pollard-machinerie n x-0 c))
      (cond ((eql 'PRIME factor) (setq factors (cons n factors) n 1))
            ((cadr factor) (setq factors (cons (setq head (car factor)) factors)
                                 n (/ n head)))
            ((null factor) (while (= (- n 2)
                                            (setq c (1+ (random (1- n)))))))
            (t (setq n (/ n (setq head (car factor)))
                     factors (append factors
                                     (rho-pollard head (random head)))))))
    factors))


(defun check-primitive (a p)
  (unless (and (miller-rabin p) (= 1 (gcd a p)))
    (return-from check-primitive))
  (let* ((phi (1- p)) (factors (sort (remove-duplicates (rho-pollard phi (random phi))) #'<)))
    (every #'(lambda (factor) (/= 1 (mod-expt a (/ phi factor) p))) factors)))
