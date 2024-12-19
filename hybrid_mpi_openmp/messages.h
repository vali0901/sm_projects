#pragma once

const int MASTER = 0;

enum SIGNAL {
    REQUEST_WORK = 0,
    SHUT_DOWN = 1,
    WORK_DONE = 2
};

enum TAGS {
    DATA = 0,
    SIGNAL = 1,
    RESULT = 2
};

class JacobiTask {
public:
    int task_id;
    struct system *sys;
};