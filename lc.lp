#theory lc {
    linear_term {
    + : 2, unary;
    - : 2, unary;
    * : 1, binary, left;
    + : 0, binary, left;
    - : 0, binary, left
    };
    assign_term {
    +  : 4, unary;
    -  : 4, unary;
    *  : 3, binary, left;
    +  : 2, binary, left;
    -  : 2, binary, left;
    .. : 1, binary, left;
    := : 0, binary, left
    };
    show_term {
    / : 0, binary, left
    };
    minimize_term {
    + : 3, unary;
    - : 3, unary;
    * : 2, binary, left;
    + : 1, binary, left;
    - : 1, binary, left;
    @ : 0, binary, left
    };

    &assign/0 : assign_term, head;
    &sum/0 : linear_term, {<=,=,>=,<,>,!=}, linear_term, body;
    &show/0 : show_term, directive;
    &distinct/0 : linear_term, body;
    &minimize/0 : minimize_term, directive
}.

