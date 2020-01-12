(ql:quickload :herwigs-cl-utilities)
(use-package 'herwig)
(ql:quickload :alexandria)
(ql:quickload :iterate)
(use-package 'iter)

(defun make-poly (&rest coefficients)
  (cond ((null coefficients) (list 0))
        ((consp (car coefficients)) (car coefficients)) ; prevent nesting
        (t coefficients)))

(defun degree (p)
  (1- (length p)))

(defun nth-coefficient (n p)
  (alexandria:if-let ((a (nth n p)))
    a
    0))

(defun +p (&rest rest)
  (let ((internal-+ (lambda (p q)
                      (iter
                        (for pk in p)
                        (for qk in q)
                        (for i upfrom 1)
                        (collect (+ pk qk) into acc)
                        (finally (return (append acc (subseq p i) (subseq q i))))))))
    (reduce internal-+ (mapcar #'mklist rest) :initial-value (make-poly 0))))

(defun -p (p &rest rest)
  (+p p (mapcar #'- (apply #'+p rest))))

(defun *p (&rest rest)
  (let ((internal-*
          (lambda (p q)
            (let ((new-degree (+ (degree p) (degree q))))
              (flet ((kth-coeff (k)
                       (iter (for l from 0 to k)
                         (sum (* (nth-coefficient l p)
                                 (nth-coefficient (- k l) q))))))
                (iter (for k from 0 to new-degree)
                  (collect (kth-coeff k))))))))
    (reduce internal-* (mapcar #'mklist rest))))

(defun deriv-poly (p)
  (iter
    (for k upfrom 1)
    (for a in (cdr p))
    (collect (* k a))))

(defun eval-poly (p x)
  (reduce (lambda (this-coeff higher-terms)
              (+ this-coeff (* x higher-terms)))
          p
          :initial-value 0
          :from-end t))
