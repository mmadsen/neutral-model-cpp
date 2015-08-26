#include <spdlog/spdlog.h>
#include "timer.h"

namespace CTModels {
std::shared_ptr<spdlog::logger> clog = spdlog::stdout_logger_mt("console");
Timer timer;
};