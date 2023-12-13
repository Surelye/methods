(defpackage #:gauss
  (:use #:cl)
  (:export #:row-echelon
           #:reduced-row-echelon
           #:row-dimension
           #:column-dimension
           #:switch-rows
           #:multiply-row
           #:add-row))


(in-package #:gauss)


(declaim (inline row-dimension column-dimension))


(defun row-dimension (a)
  "Return the number of rows of a matrix A."
  (array-dimension a 0))


(defun column-dimension (a)
  "Return the number of columns of a matrix A."
  (array-dimension a 1))


(defun switch-rows (a i j)
  "Destructively swap rows I and J of a matrix A, return A."
  (dotimes (k (column-dimension a) a)
    (psetf (aref a i k) (aref a j k)
           (aref a j k) (aref a i k))))


(defun multiply-row (a i alpha)
  "Destructively multiply the Ith row of A by alpha, return A."
  (dotimes (k (column-dimension a) a)
    (setf (aref a i k) (* (aref a i k) alpha))))


(defun add-row (a i j alpha)
  "Destructively add to Ith row of A its Jth row multiplied by alpha, return A."
  (dotimes (k (column-dimension a) a)
    (incf (aref a i k) (* (aref a j k) alpha))))


(defun eliminate-column-below (a i j)
  "Assuming that a[i,j] is nonzero, eliminate nonzero entries below a[i,j]; return a."
  (loop for k from (+ i 1) below (array-dimension a 0)
        do (add-row a k i (- (/ (aref a k j) (aref a i j))))
        finally (return a)))


(defun eliminate-column-above (a i j)
  "Assuming that a[i,j] is nonzero, destructively eliminate nonzero entries above a[i,j]; return a."
  (loop for k below i
        do (add-row a k i (- (/ (aref a k j) (aref a i j))))
        finally (return a)))


(defun find-pivot-row (a i j)
  "Return the first row number starting from i having a nonzero entry in the jth column, or nil if it does not exist."
  (loop for k from i below (row-dimension a)
        unless (zerop (aref a k j))
          do (return k)
        finally (return nil)))


(defun find-pivot-column (a i j)
  "Return the first column number starting from i having a nonzero entry in the ith row, or nil if it does not exist."
  (loop for k from j below (column-dimension a)
        unless (zerop (aref a i k))
          do (return k)
        finally (return nil)))


(defun row-echelon (a)
  "Compute the row echelon form in-place without multiplying rows by numbers; return a."
  (loop with row-dimension = (row-dimension a)
        with column-dimension = (column-dimension a)
        with current-row = 0
        with current-col = 0
        while (and (< current-row row-dimension)
                   (< current-col column-dimension))
        for pivot-row = (find-pivot-row a current-row current-col)
        do (when pivot-row
             (unless (= pivot-row current-row)
               (switch-rows a pivot-row current-row))
             (eliminate-column-below a current-row current-col)
             (incf current-row))
        do (incf current-col)
        finally (return a)))


(defun reduce-row-echelon (a)
  "Assuming that a has the row echelon form, compute the reduced row echelon form in-place without multiplying rows by numbers; return a."
  (loop for i below (row-dimension a)
        for j = (find-pivot-column a i 0) then (find-pivot-column a i j)
        while j
        unless (= 1 (aref a i j))
        do (multiply-row a i (/ 1 (aref a i j)))
        do (eliminate-column-above a i j)
        finally (return a)))


(defun reduced-row-echelon (a)
  "Compute the reduced row echelon form in-place; return a."
  (reduce-row-echelon (row-echelon a)))
