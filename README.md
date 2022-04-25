	The Java Editor program contains a text editor, email editor, table editor, and image viewer.
	
	The program also includes the Math, Number, Matrix, PublicKey, Signature, Cipher, and Convert
	classes. These classes contain all the ciphers, algorithms, protocols, and software required to
	do cryptography. The Mail, PopMail, and SendMail classes contain the software required to send
	and retrieve email.
	
	All the software is contained in one file so the source code can be executed without compiling the
	program. No other packages or modules are required to run the program or to use encrypted email.
	Software developers can copy and paste these classes into other free and open source software pro-
	grams that have compatible licenses. This means that the software can be used with a more permissive
	license but not a more restrictive or proprietary license.
	
	The Mail program can send and retrieve messages using POP mail or Post Office Protocol. The Mail
	program has a test mail feature on the help menu for the user to try the mail program without being
	connected to the internet or using a real email account. The help menu of the Mail program also has
	a document called "How to use pop mail" which explains how to use the program.
	
	Imap is not included because the protocol is more complicated to implement than POP mail. Imap al-
	lows multiple users to access an email account from different computers which is useful for some
	companies or organizations that have to reply to large numbers of emails. POP mail also allows mult-
	iple users to access an email account if none of the users deletes the new messages or only the old
	messages are deleted.
	
	Imap allows users to change the state of the messages on the server, but the POP mail protocol could
	be amended or the email servers could be upgraded to include this feature. POP mail servers could
	also be upgraded to allow multiple users to retrieve and delete emails by assigning a hash value or
	time stamp to each message so that the retrieve and delete commands could use the number assigned to
	the messages instead of the ordinal / cardinal numbers that are used to enumerate the messages.
	Otherwise if multiple users list the emails and try to delete messages using the ordinal numbers,
	the email messages on the clients' computers will not correspond to messages on the server computer
	because the messages get re-numbered every time one of the users deletes a message and signs out,
	and the wrong messages will get deleted or retrieved.
	
	The email encryption program uses a composite key that has multiple public key ciphers. The public
	key agreements are reduced modulo F8 = 2 ^ 256 + 1 and then the key agreements are xor-ed to gener-
	ate the composite secret key, session key or encryption key. Each public key agreement or cipher
	functions as a one-time pad to encrypt the other public key agreements which are also one-time pads
	or ciphers. The composite key is then used to initialize a hash function that generates another
	one-time pad for the message encryption.
	
	The private key encryption uses a hashing function that encrypts the private key and the data. No
	encryption ciphers or invertible functions are used for encryption because there is no proof that
	any invertible function that is used more than once is secure, and private keys are not supposed to
	be reused for public or private key ciphers.
	
	The public key agreement or encryption is unbreakable since every public key cipher would have to
	be broken to solve for the composite secret key. Also, the program doesn't use broken ciphers such
	as RSA or the integer Diffie-Hellman cipher which are not based on any hard math problem such as
	factorization or discrete logarithms. The software includes 16 public key ciphers (including 14
	Diffie-Merkle-Hellman ciphers and 2 Merkle-Hellman / knapsack ciphers) and 1 matrix digital sig-
	nature algorithm.
	
	The email text, file attachments, and file descriptions are each encoded in base-64, and then the
	encoded data are concatenated using newline chars (\n\n), encrypted, and re-encoded in base 64 to
	remove special characters from the encryption method such as newlines, carriage returns, and end
	of message or end of file chars. This expands the data to (4/3)^2 = 16/9 the size because base-64
	encoding maps 6 bits of data to 8-bit chars and two encodings are used to package the data. (The
	public key encryption method includes a base-64 encoding because it has to encode the cipherdata to
	attach the one-time public keys.) Other protocols may use one encoding but this would only reduce
	the expansion to 3/4 the size.
	
	The text editor and email program were written to test the public key software and to show develop-
	ers how to use and implement the public key ciphers in other programs, but anybody who knows how to
	install Java and run a java program can use the program to send and receive encrypted emails.
	
	The sender and receiver have to be using the same ciphers and protocols because the software is not
	compatible with other encryption programs. Users also may have to keep upgrading their software,
	changing their public keys, and re-encrypting their files or directories if the implementation of
	the ciphers or encryption protocols changes.
	
	
	** The current version is backward compatible with the previous versions for file encryption, but
	users should decrypt and re-encrypt their directories or folders using the new program because the
	private key cipher was modified to include permutations and future versions may not be backward com-
	patible. Also, the older versions may not accept the new public key as valid for sending encrypted
	email because two new ciphers were added, but they will still be able to decrypt messages or read
	encrypted emails.
	
	A few errors were also corrected in the software so the compiler's Xlint doesn't issue warnings
	every time the program is compiled; errors in the Number and Fourier classes were corrected; errors
	in the userpass menu item and private key encryption menu item were corrected; the file compression
	was corrected so that attached files are compressed to ~ 1/4 of their size except for files that are
	incompressible; the sign out / log out method was modified so that if the email server or mail pro-
	gram becomes unresponsive or the wifi loses its connection it will end the program in a few seconds
	so the user doesn't have to close the terminal or open the System monitor to find and terminate the
	process; a redundant encoding was removed by replacing the newlines in the encrypted and encoded da-
	ta with a base-16 separator to make it base 64; the public key ciphers were rearranged; errors in
	the readMessage method were corrected so that the messages and attached files are detached and dis-
	played correctly for encrypted and unencrypted emails; a public key padding error was corrected so
	the decryption method removes the padding / space chars appended to the message; the SavedEmails
	class was modified to sort the emails in chronological order and to view, save, or delete the at-
	tached files and edit the messages; the mail class was modified to save the message states in the
	user's mail directory by clicking on the message icons; an icon / font size error was corrected; a
	file description error was corrected that caused the delete attached files to display the file de-
	scriptions in base 64 for unencrypted emails; a passphrase substring error was corrected in the mail
	class; the sendmailframe font size was adjusted to make it the same size as the retrievemail frame
	font size and a line of code was removed from the mouse wheel listener that changed the size of the
	frame instead of the font if the control button was pressed and the mouse wheel was scrolled; an
	error in the readMessageStates method was corrected; the SavedEmails class was modified to append
	the sender's from address to the saved messages even if they have attached files; the ViewSaved-
	EmailsListener class was modified so that it creates only one instance of the SavedEmails class or
	opens only one dialog box even if the user clicks more than once on the view saved emails menu item;
	a few deprecated methods such as Frame pack() and Filechooser showDialog() were replaced even though
	the compiler doesn't issue warnings for some deprecated methods because the warnings are suppressed;
	the find class was modified so that it doesn't show the number of occurrences for an empty string;
	an error that caused the message states to reset to the unread state when a message was deleted was
	corrected; the PassphraseDialog class was rewritten to extend JDialog instead of JPanel and the code
	was modified so that the modal variable is set to false so the constructor doesn't block and the pro-
	gram can use the object returned by the constructor to set the font, color, and other variables, and
	then the modality is changed to true by the readPassphrase and readDialog methods so that the Dialog
	setVisible method blocks until the user clicks the ok button and the passphrase size and email address
	are validated; the document / file type detection was corrected so the program correctly displays html
	documents instead of trying to display them as image files which caused the dialog frame to collapse;
	
	the hyperactive class was modified to copy the url address to the clipboard so the user can copy and
	paste the address into a web browser if an email provider like yandex sends messages to clients using
	html that has hyperlinks; an icon / font size error was corrected that caused different email panels
	to have different button / icon sizes set by the readMailSettings method unless the frame was resized
	for the unselected tabs or panels; the SavedEmails variable or object was moved from the RetrieveMail-
	Frame class to the EmailPanel so that different email tabs have separate saved email frames; the state
	Changed method for the RetrieveMailFrame was modified to show and hide the saved emails frames for
	different usernames or email panels if a tab is selected or unselected; the checkDelete box method was
	modified so that checking a delete box doesn't do a read all button click which caused the screen com-
	ponents to get resized every time a box was checked or unchecked and also caused the textarea setText
	method to throw an exception if a check box was checked and unchecked; the reverse colors button was
	modified so that the button is disabled while the program is listing or reading the messages; the
	listing = true and reading = true statements were moved outside of the list and read threads so that
	they get set immediately after the user clicks the list or read button or else the color button would
	still be enabled until the list or read thread is started which caused two background colors to appear 
	simultaneously on the same list panel if two email tabs were open and the user clicked the reverse
	color button while the program was listing the messages; and the PublicKey decrypt(String, byte[])
	method was modified so that it can decrypt ciphertext using any delimiter for the prepended one-time,
	transient or ephemeral public keys such as "\n\n", "-", or the base 16 chars 0 to f.
	
	
	
	
	
	
	Instructions for running Java programs on Linux
	
	(Your computer should have at least 8 GB of memory
	if you run Java and a web browser at the same time
	or else your computer could run out of memory.)
	
	
	Downloading the Java Development Kit (JDK)
	
	To download the JDK, go to jdk.java.net
	Click on the link that says  Ready for use: JDK 18.
	
	Then click on a tar file link that says tar.gz.
	Choose the correct file for your processor architecture
	which should be x64 for Intel or Aarch for Arm processor.
	
	A dialog box appears that says read or save file.
	Click on the button that says Save File.
	
`	This should download and save the file
	openjdk-18_linux-x64_bin.tar.gz
	in the Downloads folder / directory.
	
	
	
	
	Installing the Java Development Kit (JDK) and running the Java Editor program
	
	0.  Download the file openjdk-18_linux-x64_bin.tar.gz  from the website jdk.java.net/18.
	
	1.  Drag and drop or copy and paste the Editor.java file to the Downloads folder.
	
	2.  Open a terminal and copy and paste the commands or the command line
	
	    cd; sudo mkdir -p /usr/jdk; cd; sudo cp ./Downloads/openjdk-18_linux-x64_bin.tar.gz /usr/jdk;
	    cd /usr/jdk; sudo tar zxvf openjdk-18_linux-x64_bin.tar.gz; cd;
	
	    (the -p option suppresses the error message if the directory already exists and creates the parent
	    directories as needed)
	
	3.  To run the Editor program, copy the Editor.java file to the Downloads directory and type the command
	
	    cd; /usr/jdk/jdk-18/bin/java ./Downloads/Editor.java (text, table, image, mail)
	
	If you add an argument after the file name then the program will display the text editor, table editor,
	image viewer, or email editor.
	
	(The Editor program has a table editor, html viewer, and image viewer because other editors are not
	able to display encrypted files or directories. The text editor, html viewer, and image viewer programs
	don't have to decrypt and re-encrypt the Documents and Pictures folders because they only read and de-
	crypt the file input to the program. The files on the disk remain encrypted and unmodified unless the
	user decrypts them.)
	
	
	All the commands can be concatenated into a single line using the semicolon as a delimiter.
	
	If you are running a live version of Linux, you can drag and drop the openjdk-18_linux-x64_bin.tar.gz
	file and the Editor.java file to the Downloads folder from a USB device and then copy and paste the
	single command line
	
	cd; sudo mkdir -p /usr/jdk; cd; sudo cp ./Downloads/openjdk-18_linux-x64_bin.tar.gz /usr/jdk; cd /usr/jdk;
	sudo tar zxvf openjdk-18_linux-x64_bin.tar.gz; cd; /usr/jdk/jdk-18/bin/java ./Downloads/Editor.java
	
	or for the email client
	
	cd; sudo mkdir -p /usr/jdk; cd; sudo cp ./Downloads/openjdk-18_linux-x64_bin.tar.gz /usr/jdk; cd /usr/jdk;
	sudo tar zxvf openjdk-18_linux-x64_bin.tar.gz; cd; /usr/jdk/jdk-18/bin/java ./Downloads/Editor.java mail
	
	into the terminal using the Edit -> Paste command or the popup menu.
	
	If the directory or folder name has a space character, then you have to use the back slash '\' before
	the space char to escape it. For example, a folder named My Documents would be written My\ Documents in
	the command line.
	
	
	
	It is faster to compile the program once so that the program doesn't have to be re-compiled every time.
	
	If the jdk is not installed in your computer, you first have to untar the openjdk-18 using the command
	
	 cd; sudo mkdir -p /usr/jdk; cd; sudo cp ./Downloads/openjdk-18_linux-x64_bin.tar.gz /usr/jdk;
	    cd /usr/jdk; sudo tar zxvf openjdk-18_linux-x64_bin.tar.gz; cd;
	
	To compile the Editor program, copy the Editor.java file to the Downloads folder and then copy and paste
	the command line
	
	cd; mkdir -p ./EditorClassFiles; /usr/jdk/jdk-18/bin/javac -Xlint -d ./EditorClassFiles ./Downloads/Editor.java;
	
	To run the compiled Editor or Mail program, use the command
	
	cd; /usr/jdk/jdk-18/bin/java -cp ./EditorClassFiles Editor   or
	    /usr/jdk/jdk-18/bin/java -cp /home/username/EditorClassFiles Editor
	
	
	To remove or delete the jdk directory from your computer, use the command
	
	sudo rm -r -f /usr/jdk
	
	
	
	
	
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
	\n, 3 size 1 \n, ..., or  1 size timestamp state \n, 2 size timestamp state \n, 3 size timestamp
	state \n, etcetera.
	
	This would be backward compatible with the POP mail protocol because it would only display a num-
	ber if a user changes the state of a message. Also the client could retrieve and delete messages
	using the ordinal / cardinal numbers or the time stamps. If multiple users want to retrieve and
	delete the same messages simultaneously they would have to use a newer POP mail program that re-
	trieves and deletes messages using the timestamps or message hashes.
	
	The client program stores the message hashes and message states in a file but the user has to use
	the same computer or store the mail folder / directory on a USB storage device to view the message
	states. (The program uses the hash of the from address + the number of bytes because the program
	doesn't know the hashes of the email messages from the List screen.)
	
	Until the problem of storing public keys on email servers is solved, email encryption will not be-
	come widely used. A few hundred thousand to a few million people might use encryption by copying
	and pasting keys, but hundreds of millions of people will not use email encryption by finding or
	requesting users' keys and then copying and pasting the keys into an email program. This is why
	less than a hundred thousand people and maybe even fewer than fifty thousand people in the world
	use encryption programs such as gpg.
	
	
	
	
	
	Public and private key ciphers used by the software
	
	The program uses hypercomplex and hyper-dimensional ciphers for public key agreement and a hash ci-
	pher for private key encryption. The public key agreement or secret key is hashed to generate a se-
	quence of random numbers which is used as a one-time pad. The ciphertext is computed by adding the
	one-time pad to the plaintext, and then the plaintext is recovered by subtracting the one-time pad
	from the ciphertext.
	
	The hash cipher is unbreakable because cryptographic hash functions are non-invertible. Even if the
	hash function could be inverted it wouldn't break the cipher because there are 2^768 pre-images for
	each hash value, and a cryptanalyst wouldn't know which one is the correct solution.
	
	The matrix public key ciphers are variants of the equations or functions
	
	          x1  x2           k1      k2               -x2  x1   x2           -k2  k1   k2
	  Y  =  A1  A2 ,   E  =  A1   Y  A2 ,  and  Y  =  A2   A1   A2 ,   E  =  A2   Y   A2   (mod p)
	
	which are similar to the Diffie-Merkle-Hellman cipher y = a ^ x, e = y ^ k (mod p) except that they
	use matrices or hypercomplex numbers instead of integers and they use multiple variables instead of
	a single variable. These ciphers are a generalization of the Diffie-Hellman cipher because they re-
	duce to the integer cipher y = a ^ x (mod p) if x2 = 0 and A1 is a 1x1 matrix.
	
	The integer Diffie-Hellman ciphers y = a x and y = a ^ x (mod p) can be generalized to use poly-
	nomials, vectors, matrices, cubes, tesseracts, or any n-dimensional object. Some of them can also
	use hypercomplex numbers such as quaternions or octonions.
	
	All numbers are dimensional objects. A real or complex number is a point on an axis or a plane, or
	a 0-dimensional object; an array or vector is a line or a 1-dimensional object; a matrix is a
	square or rectangular array of numbers or 2-dimensional object; a cube is 3-dimensional; and a
	tesseract is 4-dimensional.
	
	A public key Y is created by using a public parameter A to encrypt a private variable X. The sim-
	plest public key function is Y = A X. This is similar to private key encryption except that A is
	a non-invertible public parameter instead of a secret message, and X is the secret message encrypt-
	ed by the public parameter A. In private key cryptography A would be a plaintext message encrypted
	by a secret key matrix X, but in public key cryptography the private key X is the plaintext message
	encrypted by the public parameter A, and the public key Y is the ciphertext, cipherdata or cipher.
	Because the encryption key A is public, the security of public key cryptography is based entirely
	on the non-invertibility of the function instead of the secrecy of the private key.
	
	A recipient who wants to receive encrypted messages computes the static public key Y = A X. A send-
	er who wants to send an encrypted message computes the one-time public key Y = A K if the private
	variables are commutative or Y = K A if K and X are non-commutative. Then the sender and recipient
	compute the same public key agreement or secret key E = A K X  or  E = K A X  using only multipli-
	cation because each of them knows either K or X. A wiretapper would have to do an inversion to
	solve for K or X, but this is a hard math problem because A is chosen to be non-invertible.
	
	The cipher Y = A X doesn't work for integers or matrices (1 x 1 or n x n dimensional objects) be-
	cause A can be inverted to solve for X = A^-1 Y; even if A is a singular matrix the equation can
	still be solved for X. But the equation can be generalized to Y = X A X so that A is non-invertible
	and immovable because matrix multiplication is not generally commutative. Multiplication is commu-
	tative only for 0-dimensional numbers. (Also, because multiplication is non-commutative, there is
	no division operation for matrices except for integers or scalars; to divide a matrix by a matrix,
	the matrix has to be pre- or post-multiplied by the inverse of the divisor, and the divisor has to
	be an invertible or non-singular matrix.)
	
	Public key ciphers can also be generalized by using multi-dimensional multiplication instead of
	one-dimensional multiplication. For example, for 2 D multiplication, matrices can be multiplied
	from left to right and from top to bottom. Ciphers can be generalized further to use multi-dimen-
	sional algebra by using points on a plane a0 + a1 i instead of points on a line, points in a cube
	a0 + a1 i + a2 j, points in a tesseract a0 + a1 i + a2 j + a3 k (which is a quaternion), or points
	in any-dimensional space or hyperspace by defining i^2 == j^2 == k^2 == 1 and i j == k, j k == i,
	k i == j, ... Matrices of multi-dimensional points such as quaternions can also use multi-dimen-
	sional arithmetic in addition to multi-dimensional algebra.
	
	Ciphers can also be generalized by using a symmetric matrix of matrices A[][] = { { A1, A2 }
	{ A2, A3 } } as a public parameter, reducing the 2x2 block matrix to a 2x1 block matrix or public
	key vector Y[] = { A1^x1 A2^x2, A2^x1 A3^x2 } where x1, x2 are the private keys, and then reducing
	the public key vector Y[] to a 1x1 block matrix or secret key E = Y[1]^k1 Y2[2]^k2 == A1^(k1 x1)
	A2^(k1 x2 + k2 x1) A3^(k2 x2).
	
	The words block and matrix are synonymous because a matrix is a rectangular block of numbers. Be-
	fore they were called matrices, rectangular arrays of numbers were referred to as blocks. A block
	is a quantity, number, or section of things dealt with as a unit, such as a block of plaintext or
	ciphertext. (J.J. Sylvester used the term matrix in 1850 to refer to a rectangular block of numbers
	because a determinant is formed from a matrix, and a matrix is something from which something else
	originates, develops, or takes form.)
	
	Nonlinear multivariate equations are difficult or impossible to solve. If solving an equation such
	as Y = X A X were as simple as diagonalizing a matrix, doing a Fourier transform, or reducing a
	matrix to echelon and row canonical form, then math programs would have functions or methods for
	solving these equations and math books would explain how to solve them. Matrix and linear algebra
	books only explain how to solve the linear equation Y == A X or A X == B by pre-multiplying by the
	inverse of A to get X = A^-1 B. Even multivariable integer equations such as Pell's equation x^2
	- d y^2 == c or 1 are unsolvable without quantum computing, and the equations used in the public
	key class are much more difficult to solve than Pell's equation.
	
	The public key class and email program were created to use these non-invertible or one-way func-
	tions as public key ciphers. The class will eventually have tens of public key ciphers or puzzles
	copied from matrix and linear algebra books. It is unlikely that these ciphers could be broken un-
	less the Diffie-Hellman problem can be solved or the public key agreement can be computed without
	breaking the public key, inverting the function, or solving the underlying math problem.
	
	The public key class uses a composite key that includes several ciphers because there is no proof
	that any one-way function is non-invertible or that the implementation is correct and because the
	methods of cryptanalysis are secret. If cryptanalysts weren't secretive, users would know which
	public key ciphers are broken and would stop using them, and cryptographers would figure out how
	to strengthen the ciphers to resist these attacks. The only way to deal with this problem is to
	use a redundancy of ciphers based on different math problems.
	
	Composite keys are a game changer because a cryptanalyst would have to break every cipher, invert
	every function, or solve every equation in the public key class to read the encrypted messages.
	The user has an advantage since only one of the ciphers has to be secure for the encryption to be
	unbreakable. Breaking a few of the ciphers doesn't get a cryptanalyst anything because breaking a
	composite key is an all-or-nothing game.
	
	The ciphers in the public key class that have a many-to-one mapping of the private key X to the
	public key Y may be unbreakable by classical and quantum computing because the solution is ambig-
	uous and the private key X is only used once. Quantum computers are unable to solve math problems
	that have ambiguous solutions because they wouldn't know which solution to solve for. This is why
	a quantum computer can only attack the factorization problem by solving the discrete log problem.
	Even if the solution is unambiguous, it doesn't mean that a quantum computer can solve it; there
	has to be an algorithm or method for solving it, or the same private key X would have to be used
	more than once with a different public parameter A such as Y1 = A1 X and Y2 = A2 X.
	
	Encryption ciphers are not used in the software because they have a one-to-one mapping (function)
	of the plaintext to ciphertext. It doesn't make sense to use an encryption cipher that has a one-
	to-one mapping because there may be quantum (and classical) algorithms for breaking all of these
	ciphers. (A quantum algorithm already exists that can search for keys for any unknown function or
	black box in sub-exponential time by trying only the square root of the number of combinations,
	and there may be another quantum or classical algorithm that can find keys in polynomial time.)
	
	The RSA / coprime root extraction cipher c = m ^ e (mod n) where (e, phi(n)) == 1 (e and phi(n) are
	coprime) is not included or allowed in the public key class because coprime root extraction is not
	equivalent to integer factorization or based on any hard math problem. The problem with this cipher
	is that there is a one-to-one correspondence of the private keys m to the public keys c which makes
	the function invertible without factoring or unmultiplying the modulus n. This is why RSA was re-
	jected for digital signature standards and for encryption.
	
	The Rabin cipher c = m ^ e (mod n) where (e, phi(n)) != 1 (e and phi are co-composite) is equivalent
	to factorization because there is a many-to-one mapping of m to c. (Michael Rabin had thought of us-
	ing coprime root extraction as a public key cipher but he knew that it wasn't equivalent to factor-
	ization.) The Rabin cipher can use any exponent e > 1 by choosing a prime factor that has the same
	number in the totient whereas the RSA cipher can only use exponents e > 2 that are coprime with the
	totient.
	
	If the message m is a perfect square < n, then the message can be encrypted and decrypted by squar-
	ing and unsquaring m modulo n. If m is a perfect cube and phi(n) is divisible by 3, then the mes-
	sage can be encrypted and decrypted by cubing and uncubing m modulo n. The root of c = m ^ e (mod
	n) still has e ^ k solutions where k is the number of factors (or prime powers) in the modulus, but
	the recipient can extract the message by inverting e modulo phi(n)/e instead of modulo phi because
	the message is a perfect square or cube in addition to a quadratic or cubic residue modulo n.
	
	The Rabin / factorization cipher and the integer discrete log cipher are not used in the public key
	class because the factorization and integer discrete log problem are susceptible to quantum comput-
	ing. The Rabin cipher was included in the public key class to test the software for asymmetrical
	public key ciphers before the Merkle-Hellman ciphers were included because the Diffie-Hellman ci-
	phers are symmetrical which means that they use the same methods for public key generation and pub-
	lic key agreement.
	
	It doesn't make sense to use the Rabin cipher because the factorization problem is broken by quantum
	computing. The Rabin cipher can never be broken because factorization will always be harder than
	multiplication, but the key size would have to be at least 1 megabit if the running time of the algo-
	rithm is O(n^3) multi-precision multiplications or O(n^4.58) single-precision multiplications or op-
	erations. The fastest algorithm could not be faster than prime number generation which requires O(
	n^4) or O(n^3.58) operations. A 1 M bit key would only require O(1) == O(n^0) multi-precision multi-
	plications or O(n^1.58 == log2(3) == log(3)/log(2)) single-precision multiplications for encryption.
	(No key size is secure for RSA because the coprime root extraction problem is completely broken.
	This means that the function can be inverted as fast as it can be computed.)
	
	If an integer cipher is not based on the integer factorization / discrete log problem, then there is
	no need to factor the modulus or solve the discrete log problem. For example, the integer cipher y =
	a ^ x, e = y ^ k == a ^ (k x) (mod p) is broken without quantum computing because the cipher is based
	on log multiplication instead of log extraction. The integer digital signature algorithm is based on
	the discrete log problem or dlp because it uses one equation for the static signature key y = a ^ x
	(mod p), another equation for the one-time signature key r = a ^ k (mod p), and a third equation for
	the signature s = k m + x r (mod p-1) where p is the base modulus and p-1 is the exponent modulus.
	This signature algorithm is also not secure because the integer discrete log problem is just as bro-
	ken as the integer factorization problem.
	
	Elliptic curve ciphers Q = k P where the points are defined by the equation y^2 == x^3 + a x + b
	(mod p) are not included in the software because the elliptic curve discrete log function has a
	periodicity which makes it susceptible to quantum computing. In addition, the complexity of elliptic
	curves makes the ciphers vulnerable to attack without solving the ecdlp or underlying math problem
	if the parameters a, b, and p are not chosen correctly, and nobody knows how to choose the parameters
	of the curves to protect against all unknown attacks.
	
	Ciphers based on polynomial factorization and error-correcting codes also are not used or included
	in the public key class because they are not secure for any key size.
	
	The Merkle-Hellman / knapsack cipher c[] = r0 a[] + r[][] s[] (mod n), b = c[] (m[] mod p) where
	c[] is the recipient's static public key and b is the sender's one-time public key is included in
	the public key class because the cipher is quantum resistant, but the cipher is not enabled by
	default in the passphrase dialog because the key size is large. The knapsack key size is around
	6 K to 8 K bits, 1 K byte, 2 K chars, or 30 to 40 lines per key. This could be a problem for public
	keys that are printed on web pages because the keys occupy a lot of space, but it doesn't matter
	for keys that are stored on email servers.
	
	The Merkle-Hellman / knapsack cipher is important in cryptography because it is the only asymmet-
	rical public key cipher or invertible one-way function other than the Rabin cipher, and it is the
	only public key cipher that uses a private modulus. The knapsack cipher was also the world's first
	quantum-resistant public key cipher. The Merkle-Hellman cipher is unbreakable but it has to be
	implemented correctly or else it doesn't work.
	
	For example, if the public random table r[][] or the private / secret key s[] equals zero, then the
	static public key can be broken unless the cipher includes small random errors that are added to
	the static key. Moreover, if the one-time public key has a small solution set, then a cryptanalyst
	can find the sender's private key m[] by solving the knapsack problem and trying all the solutions.
	(Note that the multiplier r0 can be public because a cryptanalyst can find it anyway since one of
	the equations has to start with a[0] = 1 and r[0][] = 0; a[] has to start with 1 because a[] is a
	superincreasing sequence and the knapsack density has to equal 1.)
	
	Even if a cryptanalyst knows how to solve the general or non-superincreasing knapsack problem it
	won't break the one-time public key b = c[] m[] because the cipher was designed to resist this
	attack. The ciphers use an irregular public and private key so that there is a large number of
	solutions to the knapsack problem. The cryptanalyst wouldn't know which of these solutions is the
	correct key because the search space or solution set is too large to try all the combinations, but
	the recipient can always decrypt the sender's secret key m[] because the recipient knows the pri-
	vate modulus n and the secret key s[].
