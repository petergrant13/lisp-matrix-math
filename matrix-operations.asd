;;;; matrix-operations.asd
;;;; Defines the ASDF system for the Matrix Utilities project

(asdf:defsystem "matrix-operations"
  :description "Matrix utilities for basic linear algebra operations"
  :author "Your Name"
  :license "MIT"
  :version "0.1.0"
  :serial t  ; load files in order
  :components ((:file "matrix-operations")))

(asdf:defsystem "lisp-matrix-math/test"
  :description "Unit tests for matrix-utils"
  :depends-on ("matrix-operations")
  :serial t
  :components ((:file "test/tests")))
