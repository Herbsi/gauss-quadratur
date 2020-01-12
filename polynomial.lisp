(ql:quickload :herwigs-cl-utilities)
(use-package 'herwig)
(ql:quickload :alexandria)
(ql:quickload :iterate)
(use-package 'iter)

(defun make-poly (&rest coefficients)
  (if (null coefficients)
      (list 0)
      (mklist coefficients)))

(defun degree (p)
  (1- (length p)))

(defun nth-coefficient (n p)
  (alexandria:if-let ((a (nth n p)))
    a
    0))

(defun +p (&rest rest)
  (flet ((internal-+ (p q)
           (iter
            (for pk in p)
            (for qk in q)
            (for i upfrom 1)
            (collect (+ pk qk) into acc)
            (finally (return (append acc (subseq p i) (subseq q i)))))))
    (reduce #'internal-+ (mapcar #'make-poly rest) :initial-value (make-poly 0))))

(defun -p (p &rest rest)
  (+p p (mapcar #'- (apply #'+p rest))))

(defun *p (&rest rest)
  (flet ((internal-* (p q)
           (let ((new-degree (+ (degree p) (degree q))))
             (flet ((kth-coeff (k)
                      (iter (for l from 0 to k)
                            (sum (* (nth-coefficient l p)
                                    (nth-coefficient (- k l) q))))))
               (map0-n #'kth-coeff new-degree)))))
    (reduce #'internal-* (mapcar #'make-poly rest) :initial-value (make-poly 1))))

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
