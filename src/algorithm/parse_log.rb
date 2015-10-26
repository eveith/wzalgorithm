#!/usr/bin/env ruby

puts "epoch,x,y,z"

ARGF.each_line do |line|
  line.chomp!

  line =~ /Iteration\(epoch = (\d+), success = [-.\d]+, targetSuccess = [-.\d]+, lastSuccess = \d+, population = \((.+)\)\Z/ or next

  epoch = $1
  population = $2.split(',')

  # 1, 2, 5 (6 stk.)
  while population.size > 0
    puts "#{epoch},#{population[1][/[-\d.]+/]},#{population[2][/[-\d.]+/]},#{population[5][/[-\d.]+/]}"
    population.slice!(0, 6)
  end
end
