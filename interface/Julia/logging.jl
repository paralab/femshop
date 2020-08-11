#=
# Optional log writing
=#
export printerr, init_log, log_entry, log_dump_config, log_dump_prob, close_log

function printerr(msg)
    log_entry("Error: "*msg);
    println("Error: "*msg);
end

function init_log(name, dir)
    global log_file = dir*"\\"*name*".txt";
    global use_log = true;
    file = open(log_file, "w");
    println(file, "######################################");
    println(file, "# Femshop Log for: "*project_name);
    println(file, "######################################");
    println(file, "");
    close(file)
end

function log_entry(text)
    global log_line_index;
    global use_log;
    if use_log
        file = open(log_file, "a");
        println(file, string(log_line_index)*".\t"*text);
        log_line_index += 1;
        close(file);
    end
end

function log_dump_config(c=config)
    global log_line_index;
    global use_log;
    if use_log
        file = open(log_file, "a");
        println(file, string(log_line_index)*".\tDumping configuration:");
        log_line_index += 1;
        for f in fieldnames(Femshop_config)
            println(file, string(log_line_index)*".\t\t"*string(f)*" = "*string(getfield(c, f)));
            log_line_index += 1;
        end
        close(file);
    end
end

function log_dump_prob(p = prob)
    global log_line_index;
    global use_log;
    if use_log
        file = open(log_file, "a");
        println(file, string(log_line_index)*".\tDumping problem:");
        log_line_index += 1;
        for f in fieldnames(Femshop_prob)
            println(file, string(log_line_index)*".\t\t"*string(f)*" = "*string(getfield(p, f)));
            log_line_index += 1;
        end
        close(file);
    end
end

function close_log()
    global use_log;
    if use_log
        log_entry("Completed. Closing Log.");
        use_log = false;
    end
end