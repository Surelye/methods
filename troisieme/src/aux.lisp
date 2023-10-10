(defpackage #:aux
  (:use :cl))

(in-package #:aux)


(defmacro while (condition &body body)
  `(loop while ,condition
         do (progn ,@body)))


(defun mod-expt (base power divisor)
  (setq base (mod base divisor))
  (do ((product 1)) ((zerop power) product)
    (do () ((oddp power))
      (setq base (mod (* base base) divisor)
            power (ash power -1)))
    (setq product (mod (* product base) divisor)
          power (1- power))))


(defun compute-jacobi-machinerie (a b)
  (let ((r 1) (t-val) (c))
    (when (< a 0) (setq a (- a))
      (when (= 3 (mod b 4)) (setq r (- r))))
    (tagbody eliminate-evenness
       (setq t-val 0)
       (while (zerop (logand a 1))
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
  (when (and (integerp a) (integerp b) (= 1 (logand b 1)) (> b 2))
    (if (= 1 (gcd a b))
        (compute-jacobi-machinerie a b)
        0)))
