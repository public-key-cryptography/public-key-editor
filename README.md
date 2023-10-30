	The Java Editor program contains a text editor, email editor, table editor, and image viewer. The
	program also includes the Math, Number, Matrix, PublicKey, Signature, Cipher, and Convert classes.
	These classes contain all the ciphers, algorithms, protocols, and software required to do cryptog-
	raphy. The Mail, PopMail, and SendMail classes contain the software required to send and retrieve
	email.
	
	All the software is included in one file so the source code can be executed without compiling the
	program. No other packages or modules are required to run the program or to use encrypted email.
	Software developers can copy and paste these classes into other free and open source software pro-
	grams that have compatible licenses. This means that the software can be used with a more permissive
	license but not a more restrictive or proprietary license.
	
	The Mail program can send and retrieve messages using POP mail or Post Office Protocol. The Mail pro-
	gram has a test mail feature on the help menu for the user to try the mail program without being con-
	nected to the internet or using a real email account. The help menu of the Mail program also has a
	document called "How to use pop mail" which explains how to use the program.
	
	Imap is not included because the protocol is more complicated to implement than POP mail. Imap allows
	multiple users to access an email account from different computers which is useful for some companies
	or organizations that have to reply to large numbers of emails. POP mail also allows multiple users
	to access an email account if none of the users deletes the new messages or only the old messages are
	deleted.
	
	Imap allows users to change the state of the messages on the server, but the POP mail protocol could
	be amended or the email servers could be upgraded to include this feature. POP mail servers could
	also be upgraded to allow multiple users to retrieve and delete emails by assigning a hash value or
	time stamp in milliseconds to each message so the retrieve and delete commands could use the number
	assigned to the messages instead of the ordinal / cardinal numbers that are used to enumerate the
	messages. Otherwise if multiple users list the emails and try to delete messages using the ordinal
	numbers, the email messages on the clients' computers will not correspond to messages on the server
	computer because the messages get re-numbered every time one of the users deletes a message and signs
	out, and the wrong messages will get deleted or retrieved.
	
	The email encryption program uses a composite key that has multiple public key ciphers. The public key
	agreements are reduced modulo F8 = 2 ^ 256 + 1 and then the key agreements are xor-ed to generate the
	composite secret key, session key or encryption key. Each public key agreement or cipher functions as
	a one-time pad to encrypt the other public key agreements which are also one-time pads or ciphers. The
	composite key is then used to initialize a hash function that generates another one-time pad for the
	message encryption.
	
	The public keys are based on the Diffie-Hellman ciphers  Y = A X,  Y = X^-1 A X, and Y = x A1^x1 A0
	A2^x2 (mod p), where A, A0, A1, A2, and p are public parameters and x1, x2, and X are private keys.
	The equations use polynomials, vectors, matrices, cubes, and tesseracts. The vector cross product ci-
	pher Y = A (x) X,  E = Y * K uses a parallelogram as a public key and a parallelepiped as a shared
	secret key or key agreement. The matrix product cipher  Y = |A1|X1| uses 2-dimensional arithmetic
	which multiplies from left to right and from top to bottom. |X2|    The matrix polynomial discrete
	log cipher uses Y = A^x C B^x + ... + A^0 C B^0 (mod p). These equations were selected for public key
	cryptography because nonlinear, multivariate, multi-dimensional, modular and non-modular equations
	are unsolvable.
	
	The public keys also include the Rabin / factorization cipher c = m ^ 2 (mod n) where n is the cipher
	or static public key, m is the sender's secret key and c is the one-time public key; and the Merkle-
	Hellman / knapsack cipher c[] = a[] s0 + r[][] s[] (mod n), b = c[] (m[] + e[]) where the vector c is
	the cipher or static public key, a is a sequence of superincreasing integers, r is a public random ma-
	trix, s0, s, and n are secret keys, e is a vector of small random errors, and b is the sender's one-
	time public key.
	
	Messages are encrypted by choosing a random number or one-time encryption key (using the passphrase,
	the plaintext hash, and the system nano time as sources of entropy), hashing the random number to
	create a one-time pad, xor-ing the one-time pad and the plaindata or plaintext to generate the cipher-
	data or ciphertext, and then using the passphrase hash or shared secret key as a re-usable pad to en-
	crypt the random number or one-time encryption key. The receiver decrypts a message by xor-ing the en-
	crypted random number using the shared secret key, hashing the random number to create the one-time
	pad, and then xor-ing the one-time pad and the cipherdata to recover the plaindata.
	
	The public key agreement or encryption is unbreakable since every public key cipher would have to be
	broken to solve for the composite secret key. Also, the program doesn't use broken ciphers such as
	RSA or the integer Diffie-Hellman cipher which are not based on any hard math problem. The software
	includes 16 Diffie-Merkle-Hellman ciphers and 2 matrix digital signature algorithms.
	
	If any of these ciphers can be broken it will just get replaced. For example, if a cipher can be
	solved because it uses integers and a single equation, then it can be replaced by another cipher that
	uses matrices, polynomials, powers of a matrix (cube or tesseract), or fractional numbers instead of
	integers. Similarly, if a cipher can be broken because a matrix is diagonalizable or real and symme-
	tric, then it can be replaced by a non-diagonalizable, non-real, or non-symmetric matrix.
	
	The email text, file attachments, and file descriptions are each encoded in base-64, and then the en-
	coded data are concatenated using newline chars (\n\n), encrypted, and re-encoded in base 64 to remove
	special characters from the encryption method such as newlines, carriage returns, and end of message
	or end of file chars. This expands the data to (4/3)^2 = 16/9 the size because base-64 encoding maps
	6 bits of data to 8-bit chars and two encodings are used to package the data. (The public key encryp-
	tion method includes a base-64 encoding because it has to encode the cipherdata to attach the one-time
	public keys.) Other protocols may use one encoding but this would only reduce the expansion to 3/4 the
	size.
	
	The text editor and email program were written to test the public key software and to show developers
	how to use and implement the public key ciphers in other programs. The text editor was also used to
	write, test, and debug all the software except for the first few hundred lines which had to be written
	using a different text editor. Anybody who knows how to install Java and run a java program (either a
	java source code file or a java jar file) can use the program to do text editing or to send and re-
	ceive encrypted emails.
	
	The sender and receiver have to be using the same ciphers and protocols because the software is not
	compatible with other encryption programs. Users also may have to upgrade their software, change their
	public keys, and re-encrypt their files or directories if the implementation of the ciphers or encryp-
	tion protocols changes.
	
	
	** Note that the new version of the software will not decrypt email messages encrypted by previous
	versions of the software for keys that used more than 8 ciphers because an error was corrected in the
	secret key agreement method to reduce the secret key modulo F8 = 2^256 + 1 == 16^64 + 1; but messages
	saved on the user's computer will still be decryptable because file encryption uses private key cryp-
	tography; the format of the Merkle-Hellman / knapsack cipher was changed because the toString method
	was padding the one-time public key array with five or more zeros instead of one or more zeros,	one
	of the secret variables was redefined, the cipher was modified so that the lower elements of the
	static public key are permutated, the highest bit of the multiplier of the superincreasing sequence
	a[] was set to ensure that the product a[] s0 is reduced modulo n for all values of a[i] = 2^i; the
	isValidKey method in the public key class was redefined to only require a minimum of four valid ci-
	phers because some users may receive foreign keys from other users who have enabled large ciphers such
	as the factorization cipher or who may have a newer or different version of the software that contains
	additional ciphers; this allows users to send encrypted messages using a subset of the ciphers in the
	recipient's public key; and some encoding that was used to make the key digits look random was removed
	from three of the public keys because the methods had been modified to pack (and unpack) the 5-bit co-
	efficients into a number using adds and left shifts (and right shifts and subtracts) in the polynomial
	ciphers and the Number(Number[], int radix) constructor in the binary X A X public key (and the con-
	structor inverse toArray(int radix) in the secret key method) made the encoding unnessary because
	there are no non-random bits; in the original version of the public key class the elements were con-
	verted to string and concatenated using StringBuilder instead of shifting and adding them to a number
	and then converting the number to string.
	
	
	
	
	
	
	Instructions for running java programs on Linux
	
	(Your computer should have at least 8 GB of memory
	if you run Java and a web browser at the same time
	or else your computer could run out of memory.)
	
	
	Downloading the java development kit (JDK)
	
	To download the JDK, go to jdk.java.net
	Click on the link that says  Ready for use: JDK 21.
	
	Then click on a tar file link that says tar.gz.
	Choose the correct file for your processor architecture
	which should be x64 for Intel or Aarch for Arm processor.
	
	A dialog box appears that says read or save file.
	Click on the button that says Save File.
	
	This should download and save the file
	openjdk-21_linux-x64_bin.tar.gz
	in the Downloads folder / directory.
	
	
	
	Installing the java development kit (JDK) and running the java editor program
	
	0.  Download the file openjdk-21_linux-x64_bin.tar.gz  from the website jdk.java.net/21.
	
	1.  Drag and drop or copy and paste the Editor.java file to the Downloads folder.
	
	2.  Open a terminal and copy and paste the commands or the command line
	
	    cd; sudo mkdir -p /usr/jdk; cd; sudo cp ./Downloads/openjdk-21_linux-x64_bin.tar.gz /usr/jdk;
	    cd /usr/jdk; sudo tar zxvf openjdk-21_linux-x64_bin.tar.gz; cd;
	
	    (the -p option suppresses the error message if the directory already exists and creates the parent
	    directories as needed)
	
	3.  To run the Editor program, copy the Editor.java file to the Downloads directory and type the command
	
	    cd; /usr/jdk/jdk-21/bin/java ./Downloads/Editor.java (text, table, image, mail)
	
	If you add an argument after the file name then the program will display the text editor, table editor,
	image viewer, or email editor.
	
	(The Editor program has a table editor, html viewer, and image viewer because other editors are not
	able to display encrypted files or directories. The text editor, html viewer, and image viewer programs
	don't have to decrypt and re-encrypt the Documents and Pictures folders because they only read and de-
	crypt the file input to the program. The files on the disk remain encrypted and unmodified unless the
	user decrypts them.)
	
	
	All the commands can be concatenated into a single line using the semicolon as a delimiter.
	
	If you are running a live version of Linux, you can drag and drop the openjdk-21_linux-x64_bin.tar.gz
	file and the Editor.java file to the Downloads folder from a USB device and then copy and paste the
	single command line
	
	cd; sudo mkdir -p /usr/jdk; cd; sudo cp ./Downloads/openjdk-21_linux-x64_bin.tar.gz /usr/jdk; cd /usr/jdk;
	sudo tar zxvf openjdk-21_linux-x64_bin.tar.gz; cd; /usr/jdk/jdk-21/bin/java ./Downloads/Editor.java
	
	or for the email client
	
	cd; sudo mkdir -p /usr/jdk; cd; sudo cp ./Downloads/openjdk-21_linux-x64_bin.tar.gz /usr/jdk; cd /usr/jdk;
	sudo tar zxvf openjdk-21_linux-x64_bin.tar.gz; cd; /usr/jdk/jdk-21/bin/java ./Downloads/Editor.java mail
	
	into the terminal using the Edit -> Paste command or the popup menu.
	
	If the directory or folder name has a space character, then you have to use the back slash '\' before
	the space char to escape it. For example, a folder named My Documents would be written My\ Documents in
	the command line.
	
	
	
	Compiling the source code
	
	It is faster to compile the program once so that the program doesn't have to be re-compiled every time.
	
	If the jdk is not installed in your computer, you first have to untar the openjdk-21 using the command
	
	cd; sudo mkdir -p /usr/jdk; cd; sudo cp ./Downloads/openjdk-21_linux-x64_bin.tar.gz /usr/jdk;
	cd /usr/jdk; sudo tar zxvf openjdk-21_linux-x64_bin.tar.gz; cd;
	
	To compile the Editor program, copy the Editor.java file to the Downloads folder and then copy and paste
	the command line
	
	cd; mkdir -p ./EditorClassFiles; /usr/jdk/jdk-21/bin/javac -Xlint -d ./EditorClassFiles ./Downloads/Editor.java;
	
	To run the compiled Editor or Mail program, use the command
	
	cd; /usr/jdk/jdk-21/bin/java -cp ./EditorClassFiles Editor   or
	    /usr/jdk/jdk-21/bin/java -cp /home/username/EditorClassFiles Editor
	
	
	To remove or delete the jdk directory from your computer, use the command
	
	sudo rm -r -f /usr/jdk
	
	The --recursive option is required because the file is a directory and the
	rm command doesn't delete directories by default; the user has to confirm
	that the file to be removed is a directory by specifying the recursive
	option so users cannot inadvertently delete a directory instead of a file.
	
	The --force option tells the command not to prompt the user for a confir-
	mation before deleting each file and subdirectory, and it ignores nonexis-
	tent files and arguments which means that it will not inform the user that
	it cannot remove the file if there is no such file or directory.
	
	The sudo command is required because only the superuser can add or remove
	files that are not located in the user's home directory.
	
	
	
	
	Creating a compiled / executable java .jar file
	
	You can create a java archive or java jar file
	so the file doesn't have to be compiled each time.
	
	If the Editor.java file is in the Downloads folder, use the commands
	
	/usr/jdk/jdk-21/bin/javac -d TempDirectory Downloads/Editor.java;
	/usr/jdk/jdk-21/bin/jar cvf Editor.jar -C TempDirectory .;
	cd; echo "Main-Class: Editor" > temp.txt;
	/usr/jdk/jdk-21/bin/jar -u -f Editor.jar -m temp.txt;
	rm -r -f TempDirectory; rm temp.txt;
	
	by copying and pasting into the terminal.
	
	This command creates the jar file by compiling the program, creating a temporary directory
	for the compiled code or class files, creating a jar file and loading the class files into
	the jar file, creating a manifest file and saving the text "Main-Class: Editor", updating
	the jar file to include the manifest file, and then deleting the temporary class files
	directory and manifest file.
	
	The five command lines are printed using newline chars for readability (and to keep the
	horizontal scroll bar from expanding), but you could delete the four newline chars and
	four tab chars in the terminal and replace them with single space chars before executing
	the command to make it easier to scroll through the command history using the up and down
	arrow keys.
	
	(It doesn't matter if you run this command more than once because it just re-creates the
	jar file, but the new jar file may not have the same hash value as the previous jar file
	because it may include a time stamp.)
	
	
	
	The Editor.jar file will save around 2 to 6 seconds each time the program is executed
	depending on the speed of the processor or computer. Even on a laptop computer it only
	takes about 3 seconds to start the program except for very slow processors that could
	take up to 10 seconds.
	
	If you delete the compiled classes directory using rm -r -f EditorClassFiles, the Editor
	.jar file will still execute because the classes were loaded into the file. The class files
	are only required to create a jar file, not to run the file (unless the jar file is created
	to use the directory by specifying the class path in the manifest but then the jar file
	would only work on the computer on which the code was compiled).
	
	Note that the cd command can be omitted because it just changes the directory to the home
	directory. This is useful if the next command contains a relative path name or path that
	doesn't start with a slash /, but it is redundant to use cd if the next command has an
	absolute path name because then it doesn't do anything.
	
	The path name /usr/jdk/jdk-21/bin/java can be replaced by the file name java if the
	terminal knows where to find the java command. The path name is included because some
	users may be running a live version of Linux.
	
	
	
	Running the jar file
	
	The jar file can be run using the command
	
	cd; /usr/jdk/jdk-21/bin/java -jar Editor.jar (text) for the text editor, or
	
	cd; /usr/jdk/jdk-21/bin/java -jar Editor.jar mail
	cd; /usr/jdk/jdk-21/bin/java -jar Editor.jar table
	cd; /usr/jdk/jdk-21/bin/java -jar Editor.jar image
	
	for the email client, table editor, or image viewer, and the file will be executed
	immediately because the jar file contains the compiled classes or executable byte
	code instead of the source code. The java virtual machine will convert the byte code
	to the user's binary machine code depending on the user's computer or processor
	architecture.
	
	If you execute a .jar file instead of a .java file, then you have to remember to
	create a new jar file for each new version of the source code. You can do this either
	by copying and pasting the single command for creating the jar file each time, or by
	using the up arrow key on the keyboard to search the command history on the terminal
	until you find the command for creating the jar file and then pressing enter.
	
	
	
	
	
	
	
	Exchanging public keys and sending encrypted email
	
	To encrypt an email message the user has to find the recipient's public key, copy the public key
	to the clipboard, move the focus to the To: field (in the Send mail frame) by using the mouse or
	the tab key, and then press the enter key. If the To: field has an email address and it matches
	the address on the clipboard key, then a public key icon will appear next to the To: field to show
	that the key will be used to encrypt the message. The program can also fill in an empty To: field
	using the clipboard key address. If the clipboard key has no address, then a dialog box will ap-
	pear displaying the public key hash to confirm if the user wants to use the clipboard key as the
	recipient's public key. If the To: field already has an address and the clipboard key matches the
	address, then clicking the Send button will also display the public key icon.
	
	The problem of distributing or exchanging public keys is solvable by the email service providers
	and email server software developers. Email server programs would have to be upgraded to allow
	users to store and retrieve public keys. Users could log in to their accounts and copy and paste
	their public keys in base-16 separated by a delimiter such as '-' or the base-16 digits 0 to f,
	and then email clients could retrieve the recipient's key from the POP mail server by connecting
	to the server and then sending a request such as RETR followed by the recipient's address, or
	just sending the recipient's address to the POP mail server, and the email server would reply by
	returning the public key.
	
	For private email servers that don't have a website the server would have to allow the client /
	user to send the public key to the server after logging in to the account using a command such as
	SEND or using no command to send the key. If the client / user is logged in to a POP mail account
	and the server receives several bytes of data it would verify that the key is valid by removing
	the hyphens and testing if the data is in base 16 or only contains the chars 0 to f.
	
	Email server programs could also be upgraded so that POP mail clients could change the state of
	the messages on the server by using a POP mail command such as STAT m n where m is the message
	number and n is a state from 0 to 9. The LIST command returns an enumerated list of sizes but it
	could also return the message state number after each message size such as 1 size 0 \n, 2 size 2
	\n, 3 size 1 \n, ..., or  1 size time-stamp state \n, 2 size time-stamp state \n, 3 size time-stamp
	state \n, et cetera.
	
	This would be backward compatible with the POP mail protocol because it would only display a num-
	ber if a user changes the state of a message. Also the client could retrieve and delete messages
	using the ordinal / cardinal numbers or the time stamps. If multiple users want to retrieve and
	delete the same messages simultaneously they would have to use a newer POP mail program that re-
	trieves and deletes messages using the time stamps or message hashes.
	
	The client program stores the message hashes and message states in a file but the user has to use
	the same computer or store the mail folder / directory on a USB storage device to view the message
	states. (The program uses the hash of the from address + the number of bytes because the program
	doesn't know the hashes of the emails from the List screen or the tops of the messages.)
	
	Until the problem of storing public keys on email servers is solved, email encryption will not be-
	come widely used. A few hundred thousand to a few million people might use encryption by copying
	and pasting keys, but hundreds of millions of people will not use email encryption by finding or
	requesting users' keys and then copying and pasting the keys into an email program. This is why
	less than a hundred thousand people and maybe even fewer than fifty thousand people in the world
	use encryption programs such as gpg.
	
	
	
	
	
	Public and private key ciphers used by the software
	
	The program uses hypercomplex and hyper-dimensional ciphers for public key agreement and a hash cipher
	for private key encryption. The public key agreement or secret key is hashed to generate a sequence of
	random numbers which is used as a one-time pad. The ciphertext is computed by adding the one-time pad
	to the plaintext, and then the plaintext is recovered by subtracting the one-time pad from the cipher-
	text.
	
	The hash cipher is unbreakable because cryptographic hash functions are non-invertible. Even if the
	hash function could be inverted it wouldn't break the cipher because there are 2^768 pre-images for
	each hash value, and a cryptanalyst wouldn't know which one is the correct solution.
	
	The matrix public key ciphers are variants of the equations or functions
	
	          -x2   x1   x2          -1   x              x
	Y  =  x  A2   A1   A2 ,   Y  =  X   A   X ,  Y  =  A   X ,     Y  =  A (x) X ,   Y  =  X1  A  X2
	
	          -k2   k1   k2          -1   k              k
	E  =  k  A2   Y    A2 ,   E  =  K   Y   K ,  E  =  A   Y  K ,  E  =  Y (*) K ,   E  =  K1  Y  K2
	
	which are similar to the Diffie-Merkle-Hellman cipher y = a ^ x, e = y ^ k (mod p) except that they
	use matrices or hypercomplex numbers instead of integers and they use multiple variables instead of a
	single variable. These ciphers are a generalization of the Diffie-Hellman cipher because they reduce
	to the integer cipher y = a ^ x (mod p) if x2 = 0 and A1 is a 1x1 matrix.
	
	The integer Diffie-Hellman ciphers y = a x and y = a ^ x (mod p) can be generalized to use polynomi-
	ials, vectors, matrices, cubes, tesseracts, or any n-dimensional object. Some of them can also use
	hypercomplex numbers such as quaternions or octonions and fractional numbers or non-integers.
	
	All numbers are dimensional objects. A real or complex number is a point on an axis or a plane, or a
	0-dimensional object; an array or vector is a line or a 1-dimensional object; a matrix is a square or
	rectangular array of numbers or 2-dimensional object; a cube is 3-dimensional; and a tesseract is 4-
	dimensional.
	
	A public key Y is created by using a public parameter A to encrypt a private variable X. The simplest
	public key function is Y = A X. This is similar to private key encryption except that A is a non-
	invertible public parameter instead of a secret message, and X is the secret message encrypted by the
	public parameter A. In private key cryptography A would be a plaintext message encrypted by a secret
	key matrix X, but in public key cryptography the private key X is the plaintext message encrypted by
	the public parameter A, and the public key Y is the ciphertext, cipherdata or cipher. Because the
	encryption key A is public, the security of public key cryptography is based entirely on the non-
	invertibility of the function instead of the secrecy of the private key.
	
	A recipient who wants to receive encrypted messages computes the static public key Y = A X. This is
	the same equation as A X == B but the equation is written using X and Y instead of A and B because
	X is chosen and the equation is computed or evaluated as a function Y = f(x) = A X. If A and B are
	chosen then the equation has to be solved for X instead of evaluated for X and the equation would be
	written as A X == B instead of Y = A X. By placing the X variable on the left or right side the read-
	er knows whether the equation is being solved for X or evaluated for X.
	
	A sender who wants to send an encrypted message computes the one-time public key Z = A K if the pri-
	vate variables are commutative or Z = K A if K and X are non-commutative. Then the sender and recipi-
	ent compute the same public key agreement or secret key E = A K X  or  E = K A X  using only multi-
	plication because each of them knows either K or X. A wiretapper would have to do an inversion to
	solve for K or X, but this is a hard math problem because A is chosen to be non-invertible.
	
	The cipher Y = A X doesn't work for integers or matrices (1 x 1 or n x n dimensional objects) because
	A can be inverted to solve for X = A^-1 Y; even if A is a singular matrix the equation can still be
	solved for X. But the equation can be generalized to Y = X A X so that A is non-invertible and immov-
	able because matrix multiplication is not generally commutative. Multiplication is commutative only
	for 0-dimensional numbers such as integers, complex numbers, and quaternions which are points on a
	line, plane, or tesseract. (Also, because multiplication is non-commutative, there is no division
	operation defined for matrices except for integers or scalars; to divide a matrix by a matrix, the
	matrix has to be pre- or post-multiplied by the inverse of the divisor, and the divisor has to be an
	invertible or non-singular matrix.)
	
	Public key ciphers can also be generalized by using multi-dimensional multiplication instead of one-
	dimensional multiplication. For example, for 2-D multiplication, matrices can be multiplied from left
	to right and from top to bottom. Ciphers can be generalized further to use multi-dimensional algebra
	by using points on a plane a0 + a1 i instead of points on a line, points in a cube a0 + a1 i + a2 j,
	points in a tesseract a0 + a1 i + a2 j + a3 k, or points in any-dimensional space or hyperspace by
	defining i^2 == j^2 == k^2 == 1 and i j == k, j k == i, k i == j, ... Matrices of multi-dimensional
	points such as quaternions can also use multi-dimensional arithmetic in addition to multi-dimensional
	algebra.
	
	Ciphers can also be generalized by using a symmetric matrix of matrices such as the 2x2 block matrix
	A[][] = { { A1, A2 }, { A2, A3 } } as a public parameter, reducing the 2x2 block matrix to a 2x1 block
	matrix or public key vector Y[] = { A1^x1  A2^x2 , A2^x1  A3^x2 } where x1, x2 are the private keys,
	and then reducing the public key vector Y[] to a 1x1 block matrix or secret key E = Y[1] ^ k1 Y2[2] ^
	k2 == A1 ^ (k1 x1) A2 ^ (k1 x2 + k2 x1) A3 ^ (k2 x2).
	
	The words block and matrix are synonymous because a matrix is a rectangular block of numbers. Before
	they were called matrices, rectangular arrays of numbers were referred to as blocks. A block is a
	quantity, number, or section of things dealt with as a unit, such as a block of plaintext or cipher-
	text. (J.J. Sylvester used the term matrix in 1850 to refer to a rectangular block of numbers because
	a determinant is formed from a matrix, and a matrix is something from which something else originates,
	develops, or takes form.)
	
	Nonlinear multivariate equations are difficult or impossible to solve. If solving an equation such as
	Y = X A X were as simple as diagonalizing a matrix, doing a Fourier transform, or reducing a matrix
	to echelon and row canonical form, then math programs would have functions or methods for solving
	these equations and math books would explain how to solve them. Matrix and linear algebra books only
	explain how to solve the linear equation Y == A X or A X == B by pre-multiplying by the inverse of A
	to get X = A^-1 B. Even multivariable integer equations such as Pell's equation x^2 - d y^2 == c or 1
	are unsolvable without quantum computing, and the equations used in the public key class are much
	more difficult to solve than Pell's equation.
	
	The public key class and email program were created to use these non-invertible or one-way functions
	as public key ciphers. The class will eventually have tens of public key ciphers or puzzles copied
	from matrix and linear algebra books. It is unlikely that these ciphers could be broken unless the
	Diffie-Hellman problem can be solved or the public key agreement can be computed without breaking the
	static public key, inverting the function, or solving the underlying math problem.
	
	The public key class uses a composite key that includes several ciphers because there is no proof that
	any one-way function is non-invertible or that the implementation is correct and because the methods
	of cryptanalysis are secret. If cryptanalysts weren't secretive, users would know which public key ci-
	phers are broken and would stop using them, and cryptographers would figure out how to strengthen the
	ciphers to resist these attacks. The only way to deal with this problem is to use a redundancy of ci-
	phers based on different math problems.
	
	Composite keys are a game changer because a cryptanalyst would have to break every cipher, invert
	every function, or solve every equation in the public key class to read the encrypted messages. The
	cryptographer or user has an advantage since only one of the ciphers has to be secure for the encryp-
	tion to be unbreakable. Breaking a few of the ciphers doesn't get a cryptanalyst anything because
	breaking a composite key is an all-or-nothing game.
	
	The ciphers in the public key class that have a many-to-one mapping of the private key X to the public
	key Y may be unbreakable by classical and quantum computing because the solution is ambiguous and the
	private key X is only used once. Quantum computers are unable to solve math problems that have ambigu-
	ous solutions because they wouldn't know which solution to solve for. This is why a quantum computer
	can only attack the factorization problem by solving the discrete log problem, not by solving the dif-
	ference of squares equation x^2 == 1 (mod n). Even if the solution is unambiguous, it doesn't mean
	that a quantum computer can solve it; there has to be an algorithm or method for solving it, or the
	same private key X would have to be used more than once with a different public parameter A.
	
	Encryption ciphers are not used in the software because they have a one-to-one mapping (or function)
	of the plaintext to ciphertext. It doesn't make sense to use an encryption cipher that has a one-to-
	one mapping because there may be quantum (and classical) algorithms for breaking all of these ciphers.
	(A quantum algorithm already exists that can search for keys for any unknown function or black box in
	sub-exponential time by trying only the square root of the number of combinations, and there may be
	another quantum or classical algorithm that can find keys in polynomial time.)
	
	The RSA / coprime root extraction cipher c = m ^ e (mod n) where (e, phi(n)) == 1 (e and phi(n) are
	coprime) is not included or allowed in the public key class because coprime root extraction is not
	equivalent to integer factorization or based on any hard math problem. The problem with this cipher
	is that there is a one-to-one correspondence of the private keys m to the public keys c which makes
	the function invertible without factoring or unmultiplying the modulus n. This is why RSA was reject-
	ed for digital signature standards and for encryption. (Note that RSA refers to the cipher while co-
	prime root extraction refers to the underlying math problem on which the cipher is based.)
	
	The Rabin cipher c = m ^ e (mod n) where (e, phi(n)) != 1 (e and phi are co-composite such as e =
	2^k) is equivalent to factorization because there is a many-to-one mapping of m to c. (Michael Rabin
	had thought of using coprime root extraction as a public key cipher but he knew that it wasn't equiv-
	alent to factorization.) The Rabin cipher can use any exponent e > 1 by choosing a prime factor that
	has the same number in the totient whereas the RSA cipher can only use exponents e > 2 that are co-
	prime with the totient.
	
	If the message m is a perfect square < n, then the message can be encrypted and decrypted by squaring
	and unsquaring m modulo n. If m is a perfect cube and phi(n) is divisible by 3, then the message can
	be encrypted and decrypted by cubing and uncubing m modulo n. The root of c = m ^ e (mod n) still has
	e ^ k solutions where k is the number of factors (or prime powers) in the modulus, but the recipient
	can extract the message by inverting e modulo phi(n)/e instead of modulo phi because the message is a
	perfect square or cube in addition to a quadratic or cubic residue modulo n.
	
	The Rabin / factorization cipher is susceptible to quantum and classical computing because there are
	sub-exponential, quantum, and polynomial-time algorithms for factoring integers. For the cipher to be
	secure or unbreakable, the key size has to be on the order of 10^5 bits if the running time of the
	factorization algorithm is on the order of O(n^2) exponentiations, O(n^3) multiplications, or O(n^
	4.58) operations where n is the number of bits. The Rabin / factorization cipher is included in the
	public key class but it is not enabled by default because the key size is large. (The integer discrete
	log cipher y = a ^ x mod n is not included in the public key class because the cipher requires an ex-
	ponentiation to compute the public key y instead of a multiplication or squaring to compute the one-
	time public key c = m^2 mod n for factorization.)
	
	Elliptic curve ciphers Q = k P where the points are defined by the equation y^2 == x^3 + a x + b mod p
	are not included in the software because the elliptic curve discrete log function has a periodicity
	which makes it susceptible to quantum computing. In addition, the complexity of elliptic curves makes
	the ciphers vulnerable to attack without solving the ecdlp or underlying math problem if the parame-
	ters a, b, and p are not chosen correctly, and nobody knows how to choose the parameters of the curves
	to protect against all unknown attacks. Many people are suspicious or distrustful of elliptic curve
	ciphers because the equations are complicated and they have a large attack surface.
	
	In 2021 we wrote that elliptic curve ciphers that are based on isogenies are quantum resistant but are
	almost certainly broken since they are being approved for standardization and cryptanalysts have had
	over a decade to study them. Just because a cipher is quantum resistant doesn't mean that the cipher
	is also classical resistant or resistant to classical computing.
	
	In 2022, after we wrote that the cipher was almost certainly broken because it had been approved for
	standardization and it was being promoted and backed by a few companies, a method was published for
	breaking the supersingular isogeny key exchange cipher. If the authors hadn't published their paper,
	this algorithm would have been standardized and implemented in software programs along with the other
	broken encryption ciphers, including polynomial factorization, error-correcting code ciphers, and the
	learning with errors or LWE cipher.
	
	This example shows that the reason for the cipher competition is to discover which ciphers or equa-
	tions are complicated enough that only a few mathematicians or cryptanalysts can break or solve them,
	and then standardize those broken ciphers. This mistake or embarrassment occurred because in approving
	this cipher they underestimated the number of mathematicians who can comprehend the math that was used
	to break the cipher. For the methods of cryptanalysis to remain secret, the number of mathematicians
	who can break a cipher has to be in the single digits (such as for solving coprime root extraction,
	factorization, and the integer Diffie-Hellman problem) and in this case the number was in the double
	digits because there are tens of mathematicians who can understand the math for breaking supersingular
	isogeny key exchange.
	
	Competitions are good for many things but public key cryptography is not one of them because it just
	selects ciphers, functions, or equations that only a few people in the world know how to break, in-
	vert, or solve, and it gives users a false sense of security and confidence in the ciphers. Some users
	reassure themselves that because ciphers such as coprime root extraction or RSA have withstood many
	decades of public cryptanalysis, that this gives them a certain level of confidence in the security of
	the ciphers which is a false or erroneous assumption because cryptanalysts are secretive.
	
	Another broken cipher that is being backed by a number of companies is the learning with errors ci-
	pher. In the LWE cipher, the recipient chooses a prime (or prime power) modulus q, a public array a[],
	a private key s, and a secret random error array e[] where the sum of the elements is smaller than
	q/2, and then computes the static public key b[] = a[] s + e[] modulo q. (The random errors can be
	discarded because they are not used for decryption.) To encrypt a binary message m[], for each bit or
	element in m the sender chooses a subset of elements in a and b and then calculates u = the subset sum
	of the a elements mod q, and v = the subset sum of b[] + [q/2] m[i] (mod q) where q/2 m[i] is either 0
	or q/2. For each bit m[i] the one-time public key is the pair (u, v), and for m[] the one-time public
	key is the array of pairs or doubles { { u0, v0 }, { u1, v1 }, { u2, v2 }, ... }.
	
	Since a[] s + e[] == b[], multiplying u (== a subset sum of a[]) by s approximately equals v (== a
	subset sum of b[] + 0 or q/2); therefore the recipient can use the private / secret key s to decrypt
	or recover the message bit by calculating m == (v - s u modulo q) / [q/2] because the difference v -
	s u (mod q) equals 0 or q/2 plus the sum of the errors which is not large enough to change the quo-
	tient. But a cryptanalyst who knows how to solve the subset sum problem can also decrypt the message
	by inverting u to find the indices of the samples and then using the indices to find the subset sum of
	b[] and solving for m[i] == (v - the subset of b[] (mod q)) / [q/2]. Even if the subset sum problem
	has a many-to-one mapping, any solution to the subset sum problem will break the cipher. A cryptana-
	lyst may also be able to break the static public key because the equations are linear and the modulus
	is public unlike the knapsack cipher which is also linear and has small errors but uses a private mod-
	ulus.
	
	The Merkle-Hellman / knapsack cipher c[] = s0 a[] + r[][] s[] (mod n), b = c[] (m[] + e[]) is included
	in the public key class because it is the only cipher that uses a private modulus. Unlike the LWE ci-
	pher, this cipher is secure because it uses random errors in the static public key c[] and the one-
	time public key b. Unless the static public key could be broken, the one-time public key can never be
	broken because the solution is ambiguous and the search space or solution set is too large.
