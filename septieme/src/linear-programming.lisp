(defpackage #:lp
  (:use #:cl)
  (:export #:linear-programming))


(in-package #:lp)


(defun get-basis-decomposition (basis vec)
  (let* ((dim (length vec))
         (sle (make-array (list dim (1+ dim))))
         (slices (loop for j from 0 below dim collect (aux:array-slice basis j))))
    (do ((j 0 (1+ j))) ((= dim j) (progn
                                    (setq sle (gauss:reduced-row-echelon sle))
                                    (map 'vector #'(lambda (slice) (aref slice dim))
                                         (loop for k from 0 below dim collect (aux:array-slice sle k)))))
      (aux:change-row! sle (map 'vector #'(lambda (slice) (aref slice j)) slices) j)
      (setf (aref sle j dim) (aref vec j)))))


(defun mat-mul (A B)
  (let* ((nrowsf (gauss:row-dimension A)) (ncolss (gauss:column-dimension B))
         (nrowss (gauss:row-dimension B)) (res (make-array (list nrowsf ncolss)))
         elt)
    (do ((i 0 (1+ i))) ((= nrowsf i) res)
      (do ((j 0 (1+ j))) ((= ncolss j))
        (setq elt 0)
        (do ((k 0 (1+ k))) ((= nrowss k) (setf (aref res i j) elt))
          (setq elt (+ elt (* (aref A i k) (aref B k j)))))))))


(defun get-minor (A i j)
  (let* ((dim (gauss:row-dimension A))
         (minor (make-array (list (1- dim) (1- dim))))
         slice (iter 0) jiter)
    (do ((nrow 0 (1+ nrow))) ((= dim nrow) minor)
      (when (/= i nrow)
        (setq slice (aux:array-slice A nrow)
              jiter 0)
        (do ((ncol 0 (1+ ncol))) ((= dim ncol))
          (when (/= j ncol)
            (setf (aref minor iter jiter) (aref slice ncol))
            (setq jiter (1+ jiter))))
        (incf iter)))))


(defun det (A)
  (let ((nrows (gauss:row-dimension A)) (ncols (gauss:column-dimension A))
        (det-val 0) cur-elt)
    (when (= 1 nrows ncols)
      (return-from det (aref A 0 0)))
    (when (/= nrows ncols)
      (error "Определитель неквадратной матрицы не определён."))
    (do ((j 0 (1+ j))) ((= ncols j) det-val)
      (unless (zerop (setq cur-elt (aref A 0 j)))
        (setq det-val (+ (* cur-elt (expt -1 (mod j 2)) (det (get-minor A 0 j))) det-val))))))


(defun transpose (A)
  (let* ((nrows (gauss:row-dimension A)) (ncols (gauss:column-dimension A))
         (transposed (make-array (list ncols nrows))))
    (do ((i 0 (1+ i))) ((= nrows i) transposed)
      (do ((j 0 (1+ j))) ((= ncols j))
        (setf (aref transposed j i) (aref A i j))))))


(defun invert (A)
  (let* ((nrows (gauss:row-dimension A)) (det-val (det A))
         (inverted (make-array (list nrows nrows))))
    (do ((i 0 (1+ i))) ((= nrows i) (transpose inverted))
      (do ((j 0 (1+ j))) ((= nrows j))
        (setf (aref inverted i j) (/ (* (expt -1 (mod (+ i j) 2))
                                        (det (get-minor A i j)))
                                     det-val))))))


(defun find-solutions (eqs cs ranges)
  (let* ((state (map 'vector #'- ranges)) (len-ranges (length ranges))
         (full-ranges (map 'vector #'(lambda (range) (ash range 1)) ranges))
         (cur-symb (make-array len-ranges)) (cur-pos (1- len-ranges))
         (eq-slices (loop for j from 0 below len-ranges collect (aux:array-slice eqs j)))
         found)
    (handler-case (aux:while t
                    (when (every #'(lambda (val ineq)
                                     (<= (car ineq) val (cadr ineq)))
                                 (mapcar #'(lambda (slice)
                                             (ops:dot state slice))
                                         eq-slices) cs)
                      (setq found (cons (aux:copy-array state) found)))
                    (if (< (aref cur-symb cur-pos) (aref full-ranges cur-pos))
                        (progn
                          (incf (aref state cur-pos))
                          (incf (aref cur-symb cur-pos)))
                        (progn
                          (aux:while (= (aref full-ranges cur-pos) (aref cur-symb cur-pos))
                            (setf (aref cur-symb cur-pos) 0
                                  (aref state cur-pos) (- (aref ranges cur-pos))
                                  cur-pos (1- cur-pos)))
                          (incf (aref state cur-pos))
                          (incf (aref cur-symb cur-pos))
                          (setf cur-pos (1- len-ranges)))))
      (error ()
        (return-from find-solutions found)))))


(defun linear-programming (basis sli)
  (let* ((dim (gauss:row-dimension basis)) (reduced-basis (lll:lll basis))
         (u (loop for j from 0 below dim
                  with mul = (transpose (mat-mul reduced-basis (invert basis)))
                  collect (aux:array-slice mul j)))
         (eqs (transpose reduced-basis))
         (cs (loop for j from 0 below dim
                   for slice = (aux:array-slice sli j)
                   collect (list (aref slice 0) (aref slice (1+ dim)))))
         hs p ts R ss z-solutions)
    (destructuring-bind (ortho mus norm-sqrs) (gs:gram-schmidt reduced-basis)
      (setq hs (map 'vector #'sqrt norm-sqrs)
            p (make-array dim :initial-contents (loop for j from 0 below dim
                                                      for pair = (nth j cs)
                                                      collect (/ (+ (car pair)
                                                                    (cadr pair)) 2)))
            ts (get-basis-decomposition reduced-basis p)
            R (* (sqrt dim) (apply #'max (loop for j from 0 below dim
                                               for pair = (nth j cs)
                                               collect (/ (- (cadr pair) (car pair)) 2))))
            ss (map 'vector #'(lambda (h) (floor (/ R h))) hs)
            z-solutions (find-solutions eqs cs ss)))
    (mapcar #'(lambda (z-vec)
                (map 'vector #'(lambda (slice)
                                 (ops:dot slice z-vec))
                     u))
            z-solutions)))
