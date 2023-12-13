(defpackage #:lattice-reduction
  (:use #:cl)
  (:export #:gauss-lattice-reduction))


(in-package #:lattice-reduction)


(defun gauss-lattice-reduction (b-1 b-2)
  (let* ((b-1-copy (aux:copy-array b-1)) (b-2-copy (aux:copy-array b-2))
         (b-1-norm-sqr (ops:norm-sqr b-1-copy)) (b-2-norm-sqr (ops:norm-sqr b-2-copy))
         (b-1-2-dot (ops:dot b-1-copy b-2-copy)) a a-norm-sqr r)
    (tagbody step-1
       (setq r (floor (+ (/ b-1-2-dot b-1-norm-sqr) 1/2))
             a (ops:vec-dif b-2-copy (ops:scalar-mult r b-1-copy)))
       (when (< (setq a-norm-sqr (+ b-2-norm-sqr (* r r b-1-norm-sqr) (- (* 2 r b-1-2-dot)))) b-1-norm-sqr)
         (setq b-2-copy b-1-copy
               b-1-copy a)
         (psetq b-2-norm-sqr b-1-norm-sqr
                b-1-norm-sqr a-norm-sqr
                b-1-2-dot (- b-1-2-dot (* r b-1-norm-sqr)))
         (go step-1))
       (setq b-2-copy a)) (list b-1-copy b-2-copy)))
