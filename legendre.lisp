(ql:quickload :iterate)
(use-package 'iter)

(load "polynomial.lisp")

(defvar *legendre-polynomials* (make-hash-table))
(setf (gethash 0 *legendre-polynomials*) (make-poly 1)
      (gethash 1 *legendre-polynomials*) (make-poly 0 1))

(defun legendre-polynomial (k)
  (if (gethash k *legendre-polynomials*)
      (gethash k *legendre-polynomials*)
      (let ((new-polynomial
              (*p
               (-p (*p (make-poly 0 1)
                       (legendre-polynomial (- k 1))
                       (- (* 2 k) 1))
                   (*p (legendre-polynomial (- k 2)) (- k 1)))
               (/ 1 k))))
        (setf (gethash k *legendre-polynomials*) new-polynomial)
        new-polynomial)))

(defvar *legendre-polynomials-1st-deriv* (make-hash-table))
(defun legendre-polynomial-1st-deriv (k)
  (if (gethash k *legendre-polynomials-1st-deriv*)
      (gethash k *legendre-polynomials-1st-deriv*)
      (let ((new-polynomial (deriv-poly (gethash k *legendre-polynomials*))))
        (setf (gethash k *legendre-polynomials-1st-deriv*) new-polynomial)
        new-polynomial)))

(defun legendre (k)
  (lambda (x) (eval-poly (legendre-polynomial k) x)))

(defun legendre-1st-deriv (k)
  (lambda (x) (eval-poly (legendre-polynomial-1st-deriv k) x)))

(defun fixed-point (f x0)
  (declare (type double-float x0))
  (declare (optimize (speed 3)))
  (iter
    (for k from 1 to 10000)
    (for y0 initially x0 then yk)
    (for yk next (funcall f y0))
    (finally (return yk))))

(defun newton-method (f df x0)
  (declare (type double-float x0))
  (flet ((newton-transform (g dg)
           (lambda (x)
             (declare (type double-float x))
             (- x (/ (funcall g x) (funcall dg x))))))
    (fixed-point (newton-transform f df) x0)))

(defun legendre-roots (n)
  (let ((lp (legendre (1+ n)))
        (lp1 (legendre-1st-deriv (1+ n))))
    (iter
      (for k from 0 to n)
      (collect
          (let ((x0 (cos (/
                          (* (+ (* 4 k) 3) pi)
                          (+ (* 4 n) 6)))))
            (newton-method lp lp1 x0))))))

(defun integration-weight (roots)
  (let* ((n (- (length roots) 1))
         (lp (legendre n)))
    (iter (for xk in roots)
      (collect (/
                (* 2 (- 1 (expt xk 2)))
                (* (expt (+ n 1) 2) (expt (funcall lp xk) 2)))))))

(defun gauss-quadratur (n)
  (let* ((roots (legendre-roots n))
         (weights (integration-weight roots)))
    (lambda (f)
      (iter
        (for root in roots)
        (for weight in weights)
        (sum (* (funcall f root) weight))))))

(defconstant +actual-value+ (the double-float (- (log 27.0d0) 2)))

;; Somehow during polynomial evaluation rounding errors kreep up
;; in the last 4 decimal digits or so, preventing precision past 10e-12 unfortunately
;; I have no idea why that happens or how to prevent it

(iter
  (for n in (list 2 4 8 16))
  (for int next (funcall (gauss-quadratur n) (lambda (x)
                                               (declare (type double-float x))
                                               (log (+ x 2)))))
  (for err next (abs (- int +actual-value+)))
  (format t "n: ~2,,d     int: ~,17E     err: ~,17E~%" n int err))
