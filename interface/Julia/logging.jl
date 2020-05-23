#=
# Optional log writing
=#
export printerr, init_log, log_entry, log_dump_config, close_log

log_line_index = 1;
use_log = false;
log_file = nothing;

function printerr(msg)
    log_entry("Error: "*msg);
    println("Error: "*msg);
end

function init_log(name, dir)
    file = open(dir*"\\"*name*"_log.txt", "w");
    println(file, "######################################");
    println(file, "# Femshop Log for: "*name);
    println(file, "######################################");
    println(file, "");
    global log_file = file;
    global use_log = true;
end

function log_entry(text)
    global log_line_index;
    global use_log;
    if use_log
        println(log_file, string(log_line_index)*".\t"*text);
        log_line_index += 1;
    end
end

function log_dump_config(config)
    global log_line_index;
    global use_log;
    if use_log
        println(log_file, string(log_line_index)*".\tDumping Femshop configuration:");
        log_line_index += 1;
        for f in fieldnames(Femshop_config)
            println(log_file, string(log_line_index)*".\t\t"*string(f)*" = "*string(getfield(config, f)));
            log_line_index += 1;
        end
    end
end

function close_log()
    global use_log;
    if use_log
        log_entry("Completed. Closing Log.");
        close(log_file);
        use_log = false;
    end
end