#=
# Optional log writing
=#
export printerr, init_log, log_entry, log_dump_config, log_dump_prob, close_log

function printerr(msg)
    log_entry("Error: "*msg);
    println("Error: "*msg);
end

function init_log(name, dir)
    file = open(dir*"\\"*name*".txt", "w");
    println(file, "######################################");
    println(file, "# Femshop Log for: "*project_name);
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

function log_dump_prob(problem)
    global log_line_index;
    global use_log;
    if use_log
        println(log_file, string(log_line_index)*".\tDumping problem:");
        log_line_index += 1;
        for f in fieldnames(Femshop_prob)
            println(log_file, string(log_line_index)*".\t\t"*string(f)*" = "*string(getfield(problem, f)));
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