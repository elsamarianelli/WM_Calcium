% Helper function to map characters to numbers
function num = map_char_to_number(char)
    switch char
        case 'A'
            num = 1;
        case 'B'
            num = 2;
        case 'C'
            num = 3;
        otherwise
            error('Invalid character');
    end
end