(defpackage #:aux
  (:use #:cl)
  (:export #:while
           #:array-slice
           #:copy-array
           #:is-linear-dependent?
           #:change-row!
           #:switch-rows!))


(in-package #:aux)


(defmacro while (condition &body body)
  `(loop while ,condition
         do (progn ,@body)))


(defun array-slice (arr row)
  (let ((arr-dim (gauss:column-dimension arr)))
    (make-array arr-dim
                :displaced-to arr
                :displaced-index-offset (* row arr-dim))))


(defun copy-array (array &key
                   (element-type (array-element-type array))
                   (fill-pointer (and (array-has-fill-pointer-p array)
                                      (fill-pointer array)))
                   (adjustable (adjustable-array-p array)))
  (let ((dims (array-dimensions array)))
    (adjust-array
     (make-array dims
                 :element-type element-type :fill-pointer fill-pointer
                 :adjustable adjustable :displaced-to array)
     dims)))


(defun is-linear-dependent? (arr)
  (let ((arr-copy (copy-array arr)))
    (gauss:reduced-row-echelon arr-copy)
    (every #'zerop
           (array-slice arr-copy (1- (gauss:row-dimension arr-copy))))))


(defun change-row! (arr sub-row row)
  (let ((len-row (length sub-row)))
    (do ((j 0 (1+ j))) ((= j len-row) arr)
      (setf (aref arr row j) (aref sub-row j)))))


(defun switch-rows! (arr i j)
  (let ((num-cols (gauss:column-dimension arr)))
    (do ((k 0 (1+ k))) ((= num-cols k) arr)
      (psetf (aref arr i k) (aref arr j k)
             (aref arr j k) (aref arr i k)))))
