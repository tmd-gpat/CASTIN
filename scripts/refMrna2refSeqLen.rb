current = ""
length = 0
while l = gets do
    if l[0] == '>' then
        if current.size > 0 then
            print current
            print "\t"
            puts length
            length = 0
        end
        name = l.split(' ')[0]
        current = name[1...name.size]
    else
        length += l.size-1
    end
end
print current
print "\t"
puts length
