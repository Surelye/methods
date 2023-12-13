(defpackage #:gs
  (:use #:cl)
  (:export #:gram-schmidt))


(in-package #:gs)


(defun gram-schmidt (vecs)
  (let* ((row-dim (gauss:row-dimension vecs))
         (ortho-vecs (make-array (list row-dim (gauss:column-dimension vecs))))
         (norm-sqrs (make-array row-dim))
         (mus (make-array (list row-dim row-dim))) temp-row mu)
    (aux:change-row! ortho-vecs (setq temp-row (aux:array-slice vecs 0)) 0)
    (setf (aref norm-sqrs 0) (ops:norm-sqr temp-row))
    (do ((i 1 (1+ i))) ((= row-dim i) (list ortho-vecs mus norm-sqrs))
      (aux:change-row! ortho-vecs (setq temp-row (aux:array-slice vecs i)) i)
      (do ((j 0 (1+ j))) ((= i j) (setf (aref norm-sqrs i)
                                        (ops:norm-sqr (aux:array-slice ortho-vecs i))))
        (setq mu (/ (ops:dot temp-row (aux:array-slice ortho-vecs j)) (aref norm-sqrs j)))
        (setf (aref mus i j) mu)
        (aux:change-row! ortho-vecs
                         (ops:vec-dif (aux:array-slice ortho-vecs i)
                                      (ops:scalar-mult mu (aux:array-slice ortho-vecs j)))
                         i)))))
