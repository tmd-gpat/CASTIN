min_len = ARGV[0].to_i

if not (min_len > 0) then
  STDERR.puts "specify valid minimum length"
  exit
end

while true do
  lines = []
  for i in 0...4 do
    exit if (l = STDIN.gets).nil?
    lines << l
  end
  if lines[1].strip.length >= min_len
    for i in 0...4 do
      print lines[i]
    end
  end
end

