while l = gets do
    if l[0] == '>' then
        name = l.split(' ')[0]
        puts name[1..name.size-1]
    end
end
