% steps used in solving the problem:
%
% 1. encode the matrix vhod into a row vector using run-length encoding
% make sure the encoding always stars with 0
% the result of this step should be a row vector of lengths of strings.
% Since we are always starting with 0 and there are only two possible values we canonicalize_file_name
% drop the symbol specifier
%
% 2. Find the lengths of the code substitutions for symbols of the vector we got from step 1.
% Use the Huffman encoding algorithm. Make sure to get two vectors of lengths - the lengths
% for white pixels and lengths for black pixels.
%
% 3. Use the code substitution lengths to encode the symbols using the canonical huffman code
%
% 4. Substitute symbols given by the rle algorithm with the code substitutions
% return coded image, the compression  ratio, the vector of code substitution lengths of
% white pixels and the vector of code substitution lengths of black pixels.
 
function [izhod, R, kodBela, kodCrna] = compress(vhod)
	vhod = 1 - vhod;
	% computing the encoded rows, white symbols and their frequencies and black symbols and their frequencies
	[rl_encoded, rl_white_u, rl_black_u] = rle(vhod);
	
	% computing code substitution lengths for white and black symbols
	white_lengths = get_huff_lengths(rl_white_u);
	black_lengths = get_huff_lengths(rl_black_u);
	
	% computing canonical Huffman code substitutions for white and black symbols
	white_canonical = to_canonical(white_lengths);
  black_canonical = to_canonical(black_lengths);
	
	% encode the input
	izhod = encode(rl_encoded, white_canonical, black_canonical);
	% compute R
	R = numel(izhod)/numel(vhod);
	
	% also return code substitution lengths
	kodBela = sortrows(white_lengths, [2, 1]);
	kodCrna = sortrows(black_lengths, [2, 1]);
	
endfunction

% function [rl_encoded, rl_white, rl_black] = rle(data)
% computes and returns the runlength encoding of the data matrix_type
% the encoding of a row always starts by encoding 0 - the element specifier can be omitted
%
% the input is a row vector of bits. 1s represent black pixels, 0s represent white pixels

% the ouput are matrices with unique symbols in the first column and their frequencies in the second column -
% one for symbols representing white pixels and one for symbols representing black pixels

% and a vector of structures representing rl encoded rows

function [rl_encoded, rl_white_u, rl_black_u] = rle(data)
  % initialize matrices for the rl encoded image and white/black symbols that appear in the row
  rl_encoded = [];
  white_symbols = [];
  black_symbols = [];
  % go over rows in the data matrix
  for data_row = data'
    data_row = data_row';
    % compute the rle of the passed data
    T = find(diff([data_row(1) - 1, data_row]));
    encoded_row = [diff([T, numel(data_row) + 1])];
    % if data does not start with 0, append
    if data_row(1, 1) != 0
      encoded_row = [0, encoded_row];
    endif
    white_symbols = [white_symbols, encoded_row(1:2:end)];
    black_symbols = [black_symbols, encoded_row(2:2:end)];
    row.val = encoded_row;
    rl_encoded = [rl_encoded; row];
  endfor
  
  % computing unique symbols encoding black and white pixels
  rl_white_u = unique(white_symbols);
  rl_black_u = unique(black_symbols);

  % computing frequencies and merging with the list of unique symbols 
  white_freq = histc(white_symbols, rl_white_u);
  rl_white_u = [rl_white_u', white_freq'];
  
  black_freq = histc(black_symbols, rl_black_u);
  rl_black_u = [rl_black_u', black_freq'];
  
endfunction

% function [lengths] = get_huff_lengths(symbols)
% function for computing the Huffman encoding lengths of passed symbols with their frequencies

% the function accepts a matrix where the first column contains the symbols and the second column their frequencies
% the function returns a matrix where the first column contains the symbols and the second column contains their Huffman lengths

% this function should be called twice - once for symbols representing lengths of white pixels
% and once for symbols representing black pixels
function [lengths] = get_huff_lengths(symbols)
	
	% sort by frequencies (descending) and then by symbols (ascending)
	symbols = sortrows(symbols, [-2 1]);
	
  % compute total frequency
	freq = sum(symbols(:, 2));
	
  % create "blank" matrix in which huffman length generation will be performed
  working_mat = zeros(length(symbols) + length(symbols) - 1, 5);
  
  % hold index of next vacant row
  index_vacant = length(symbols) + 1;
  
  % add symbols and frequencies to working matrix
  working_mat(1:length(symbols), 1:2) = symbols;
  
  % mask off zero values in the frequencies column
  mask = (working_mat(:, 2) == 0);
  working_mat(mask, 2) = NaN;
  
  % initialize sum_freq (holds current sum of frequencies of symbols that are being merged)
  sum_freq = 0;
  
  while sum_freq < freq
    
    % find sign with smallest frequency (>0) and save value and index
    min_freq1 = min(working_mat(:, 2));
    
		% handling duplicate frequencies ////////////////////////////////
		
		% get vector of indices with minimum frequency
		indices_pool = find(working_mat(:, 2) == min_freq1);
		
		% reduce this set to indices of symbols composed of a minimum number of elements
		indices_pool = get_simplest(working_mat, indices_pool, 1:length(symbols));
		
		I1 = get_minimum(working_mat, indices_pool);
		% ///////////////////////////////////////////////////////////////
		
    % put sign in next vacant row
    working_mat(index_vacant, 4) = I1;
    
    % mark symbol as used
    working_mat(I1, 2) = NaN;
    % finding second symbol /////////////////////////////////////////
    
    % find sign with smallest frequency (>0) and save value and index
    min_freq2 = min(working_mat(:, 2));
    
		% handling duplicate frequencies ////////////////////////////////
		
		% get vector of indices with minimum frequency
		indices_pool = find(working_mat(:, 2) == min_freq2);
		
		% reduce this set to indices of symbols composed of a minimum number of elements
		indices_pool = get_simplest(working_mat, indices_pool, 1:length(symbols));
		
		I2 = get_minimum(working_mat, indices_pool);
		% ///////////////////////////////////////////////////////////////
		
    % put sign in next vacant row
    working_mat(index_vacant, 5) = I2;
    
    % mark symbol as used
    working_mat(I2, 2) = NaN;
    
    sum_freq = min_freq1 + min_freq2;
    
    % add frequencies and put resut into working_mat
    working_mat(index_vacant, 2) = sum_freq;
    
    % increment index of next vacant rows
    index_vacant++;
    
    % continue until frequency adds up to freq
  endwhile
  
  % build code substitution lengths //////////////////////////
  
  % go over the code lengths column from bottom up
  for  i = size(working_mat)(1):-1:length(symbols) + 1
    % get indices of composing symbols
    indices_comp = working_mat(i, 4:5);
    if any(indices_comp)
			working_mat(indices_comp, 3) = working_mat(indices_comp, 3) + working_mat(i, 3) + 1;
		endif
		
  endfor
  
	% construct the resulting matrix lengths consisting of [symbols, code substituion lengths]
  lengths = working_mat(working_mat(:, 3) > 0, 3);
	lengths = lengths(1:length(symbols),:);
	lengths = [symbols(:, 1), lengths];
  
endfunction

% function [I] = get_simplest(working_mat, indices_pool, simple_indices)
% auxilliary function that selects the symbols composed of least other symbols from symbols indexed by indices_pool
function [I] = get_simplest(working_mat, indices_pool, simple_indices)
  % find complexity of symbols indexed by indices_pool
	pool = arrayfun(@(i) get_complexity(working_mat, i, simple_indices), indices_pool);
  % join vector of complexities with vector of indices
  pool = [indices_pool, pool];
  % save indices indexing symbols with minimum comlexity
	I = pool(pool(:, 2) == min(pool(:, 2)));
endfunction

% function [res] = get_complexity(working_mat, row, simple_indices)
% auxilliary function for getting the complexity of a symbol (how many elementary symbol it is composed of)
function [res] = get_complexity(working_mat, row, simple_indices)
  % if symbol at row row is a member of indices indexing one-symbol symbols then complexity is 1
	if ismember(row, simple_indices)
		res = 1;
	else
    % else recursively get complexity of elements that the symbol is composed of
		res = get_complexity(working_mat, working_mat(row, 4), simple_indices) + get_complexity(working_mat, working_mat(row, 5), simple_indices);
	endif
endfunction

% function [res] = get_minimum(working_mat, indices_pool)
% auxilliary function for getting the symbol indexed by indices_pool that has the minimum value
function [res] = get_minimum(working_mat, indices_pool, symbols_length)
  % set all zero values in the working matrix to NaN
	working_mat(working_mat == 0) = NaN;
	% recover any 0 symbol that was masked to NaN
	working_mat(isnan(working_mat(:, 1))) = 0;
  % get minimum value indexed by any indices in the indices_pool
  % in case of composed symbols use the first symbol as reference
	minimum = min(min([working_mat(indices_pool, 4); working_mat(indices_pool, 1)]));
  % find which index indexes this minimum valued symbol
	res = -1;
	for i = indices_pool'
		if working_mat(i, 4) == minimum || working_mat(i, 1) == minimum
			res = i;
      return
		endif
	endfor
endfunction



% function [canonical_code] = to_canonical(code_subs)
% computes and returns the canonical representation of the huffman code
% lengths given in the arguments

% ////////////////////////////////////////////////////////////////////////////////////////
% the first column of the i-th row of the input matrix sorted
% should contain the symbol to which the data in the following columns applies

% the second column  of the i-th row of the input matrix sorted
% should contain the length of the code substitution for the i-th symbol (sorted)

% the input matrix sorted should be sorted first by substitution length and
% then by the symbol value (both ascending);

% this function produces the huff_white and huff_black matrices that are used as
% input for function encode.
% ////////////////////////////////////////////////////////////////////////////////////////

% this function should be called twice - 
% once for white symbols representing white pixels and then for symbols representing black pixels

function [canonical_code] = to_canonical(sorted)  
	
	% sort the matrix sorted
	sorted = sortrows(sorted, [2, 1]);
	
  % create matrix for holding the canonical huffman code
  % the first column contains the length of the code substitution
  temp = zeros(size(sorted)(1), sorted(end, 2) + 2);
  temp(1:end, 2) = sorted(:, 2);
  temp(1:end, 1) = sorted(:, 1);
  sorted = temp;
  
  % let next element always be one binary number higher
  % if binary representation shorter than original append 0s
  for i = 1:size(sorted)(1)
    if i == 1
      sorted(1, 3:end) = 0;
      continue;
    end
    % processing the next code substitution ///////////////////
    % convert vector of binary digits to string
    bin_rep = num2str(sorted(i - 1, end - sorted(i - 1, 2) + 1:end));
    bin_rep = get_next_bin(bin_rep);
    % /////////////////////////////////////////////////////////
    % Check if new code substitution is long enough and increment
    % if necessary.
    while length(bin_rep) < sorted(i, 2)
      bin_rep = [bin_rep, 0];
    end
    % inject current code substitution back into array
    sorted(i, end - length(bin_rep) + 1:end) = bin_rep; 
  endfor
  % let canonical_code be the sorted matrix
  canonical_code = sorted;
endfunction

% function [bin_rep] = get_next_bin(bin_str)
% auxilliary function for incrementing the string representation of a binary number bin_str by 1
function [bin_rep] = get_next_bin(bin_str)
  % remove whitespace
  bin_str(isspace(bin_str)) = '';
	len = length(bin_str);
  % convert to decimal value
  bin_str = bin2dec(bin_str);
  % increment
  bin_str++;
  % convert to binary string
  bin_str = dec2bin(bin_str, len);
  bin_rep = bin_str - '0';
endfunction


% function [encoded_final] = encode(rl_encoded_rows, huff_white, huff_black)
%
% ////////////////////////////////////////////////////////////////////////////////////////
% function accepts a row vector representing the run-length encoded input and two mapping matrices.
% The rows are given as structures contained in the row vector. Each structure contains a member val that contains the next row

% the first column  of the encoding matrix contains the symbol, the second column contains the length of the code substitution
% and remaining columns represent the code substitution (the number is specified by the second column)
% ////////////////////////////////////////////////////////////////////////////////////////
function [encoded_final] = encode(rl_encoded_rows, huff_white, huff_black)
	encoded_final = [];
  for row = 1:numel(rl_encoded_rows)
    rl_encoded = rl_encoded_rows(row).val;
    % separate rl_encoded into odd and even parts
    rl_encoded_white = rl_encoded(1:2:end);
    rl_encoded_black = rl_encoded(2:2:end);
    
    
    % sort huff_white by first column
    % huff_white = sortrows(huff_white, 1);
    
    % get indices of uniqe symbols as they appear in the rl_encoded_white vector
    % [~, ~, uniq] = unique(rl_encoded_white);
    [~, uniq] = ismember(rl_encoded_white', huff_white(:, 1));
    % the vector of unique elements is the same as the first column of huff_white
    % substitute the symbols with corresponding code substitutions
    % uniq contains the indices of elements in the huff_white arrayfun
    % map these indices to value is huff_white
    encoded_white = arrayfun(@(val) apply_sub(val, huff_white), uniq, "UniformOutput", false);
    
    % do the same for rl_encoded_black
    % huff_black = sortrows(huff_black, 1);
    % [~, ~, uniq] = unique(rl_encoded_black);
    
    [~, uniq] = ismember(rl_encoded_black', huff_black(:, 1));
    encoded_black = arrayfun(@(val) apply_sub(val, huff_black), uniq, "UniformOutput", false);

    % joining the vectors together
    temp = cell(length(encoded_white) + length(encoded_black), 1);
    temp(1:2:end) = encoded_white;
    temp(2:2:end) = encoded_black;
    encoded_final_temp = cell2mat(temp');
    encoded_final_temp = encoded_final_temp - '0';
    encoded_final = [encoded_final, encoded_final_temp];
  endfor
  
endfunction

% function [mapping] = apply_sub(val,huff)
% apply sub: auxilliary function to be used for mapping indices to values
function [mapping] = apply_sub(val,huff)
  % get lenght of code substitutiton
  len = huff(val, 2);
  % compute the string representation of the code substitution
  mapping = huff(val, end - len + 1:end);
  mapping = num2str(mapping);
  mapping(isspace(mapping)) = '';
endfunction