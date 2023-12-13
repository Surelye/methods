(defpackage #:ops
  (:use #:cl)
  (:export #:norm-sqr
           #:norm
           #:dot
           #:vec-dif
           #:scalar-mult))


(in-package #:ops)


(defun norm-sqr (vec)
  (reduce #'+ (map 'vector #'(lambda (b-i) (* b-i b-i)) vec)))


(defun norm (vec)
  (sqrt (norm-sqr vec)))


(defun dot (f s)
  (reduce #'+ (map 'vector #'* f s)))


(defun vec-dif (f s)
  (map 'vector #'- f s))


(defun scalar-mult (scalar vec)
  (map 'vector #'(lambda (b-i) (* scalar b-i)) vec))
